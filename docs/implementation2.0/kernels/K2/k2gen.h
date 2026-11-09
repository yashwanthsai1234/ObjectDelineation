// k2gen.h -- deterministic input generators for k2test and k2bench. Own RNG (splitmix64) and integer-only geometry
// (no <random> distributions, no libm), so a seed gives the same inputs on macOS/libc++ and Linux/libstdc++ and the
// CPU and GPU builds can compare output digests.
#pragma once

#include <algorithm>
#include <cstring>
#include <utility>
#include <vector>

#include "k2.h"
#include "k2ref.h"

namespace k2gen {

using k2::i32;
using k2::i64;
using k2::u32;
using k2::u64;
using Pt = std::pair<i64, i64>;

struct Rng {
  u64 s;
  explicit Rng(u64 seed) : s(seed) {}
  u64 next() {
    u64 z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  i64 range(i64 lo, i64 hi) { return lo + i64(next() % u64(hi - lo + 1)); }  // inclusive
  bool chance(int pct) { return i64(next() % 100) < pct; }
};

// FNV-1a over raw bytes: the digest the CPU and GPU builds must agree on.
struct Digest {
  u64 h = 0xcbf29ce484222325ull;
  void bytes(const void* p, size_t n) {
    const unsigned char* c = static_cast<const unsigned char*>(p);
    for (size_t i = 0; i < n; ++i) h = (h ^ c[i]) * 0x100000001b3ull;
  }
  template <class T>
  void vec(const std::vector<T>& v) {
    const u64 n = v.size();
    bytes(&n, sizeof n);
    if (!v.empty()) bytes(v.data(), sizeof(T) * v.size());
  }
};

// ---------------------------------------------------------------- polygons
struct Polys {
  std::vector<i32> xy, roles, W, H;
  std::vector<i64> ring_start{0}, shape_start{0};
  i64 G = 1000;
  i64 rings_in_shape = 0;
  void ring(const std::vector<Pt>& p, i32 role) {
    for (const Pt& q : p) xy.push_back(i32(q.first)), xy.push_back(i32(q.second));
    ring_start.push_back(i64(xy.size() / 2));
    roles.push_back(role);
    ++rings_in_shape;
  }
  void end_shape(i32 w = 64, i32 h = 64) {
    shape_start.push_back(i64(roles.size()));
    W.push_back(w), H.push_back(h);
    rings_in_shape = 0;
  }
  i64 n() const { return i64(W.size()); }
  gpucheck::Polygons view() const {
    gpucheck::Polygons p;
    p.n = n(), p.xy = xy.data(), p.ring_start = ring_start.data(), p.shape_start = shape_start.data();
    p.roles = roles.data(), p.W = W.data(), p.H = H.data(), p.G = G;
    return p;
  }
};

inline std::vector<Pt> random_ring(Rng& r, int n, i64 x0, i64 y0, i64 x1, i64 y1) {
  std::vector<Pt> p;
  for (int i = 0; i < n; ++i) p.push_back({r.range(x0, x1), r.range(y0, y1)});
  return p;
}

// Random points sorted by exact angle around (cx, cy) (integer half-plane + cross product, ties by distance, then by
// coordinates): a star-shaped ring, simple unless points repeat or line up with the centre. ccw flips the order.
inline std::vector<Pt> star_ring(Rng& r, int n, i64 cx, i64 cy, i64 rmin, i64 rmax, bool ccw) {
  std::vector<Pt> p;
  while (int(p.size()) < n) {
    const i64 dx = r.range(-rmax, rmax), dy = r.range(-rmax, rmax), d2 = dx * dx + dy * dy;
    if (d2 < rmin * rmin || d2 > rmax * rmax) continue;
    p.push_back({cx + dx, cy + dy});
  }
  auto half = [&](const Pt& a) { const i64 dx = a.first - cx, dy = a.second - cy; return (dy < 0 || (dy == 0 && dx < 0)) ? 1 : 0; };
  std::sort(p.begin(), p.end(), [&](const Pt& a, const Pt& b) {
    const int ha = half(a), hb = half(b);
    if (ha != hb) return ha < hb;
    const i64 ax = a.first - cx, ay = a.second - cy, bx = b.first - cx, by = b.second - cy, c = ax * by - ay * bx;
    if (c != 0) return c > 0;
    const i64 da = ax * ax + ay * ay, db = bx * bx + by * by;
    if (da != db) return da < db;
    return a < b;
  });
  if (ccw) std::reverse(p.begin(), p.end());
  return p;
}

// Axis-aligned rectangle with `sub` points per side (collinear points along the sides), clockwise on screen or not.
inline std::vector<Pt> rect_ring(i64 x0, i64 y0, i64 x1, i64 y1, int sub, bool cw) {
  std::vector<Pt> p;
  for (int i = 0; i < sub; ++i) p.push_back({x0 + (x1 - x0) * i / sub, y0});
  for (int i = 0; i < sub; ++i) p.push_back({x1, y0 + (y1 - y0) * i / sub});
  for (int i = 0; i < sub; ++i) p.push_back({x1 - (x1 - x0) * i / sub, y1});
  for (int i = 0; i < sub; ++i) p.push_back({x0, y1 - (y1 - y0) * i / sub});
  if (!cw) std::reverse(p.begin(), p.end());
  return p;
}

// One shape of a random family. Caps: <= 10 rings, <= 1000 vertices per ring.
inline void random_shape(Rng& r, Polys& P, i64 G, int max_verts) {
  const int fam = int(r.range(0, 9));
  const i64 c = G / 2;
  auto nv = [&]() { return int(r.range(3, max_verts)); };
  switch (fam) {
    case 0:  // one star ring
      P.ring(star_ring(r, nv(), c, c, G / 8, G / 2 - 1, r.chance(50)), -1);
      break;
    case 1: {  // outer star with small star holes (may overlap, nest, cross)
      P.ring(star_ring(r, nv(), c, c, G / 4, G / 2 - 1, false), -1);
      const int holes = int(r.range(1, 9));
      for (int h = 0; h < holes; ++h) {
        const i64 hx = r.range(G / 4, 3 * G / 4), hy = r.range(G / 4, 3 * G / 4);
        P.ring(star_ring(r, int(r.range(3, max_verts / 4 + 3)), hx, hy, 1, G / 10, true), 0);
      }
      break;
    }
    case 2: {  // several outers, random places (may overlap or nest)
      const int k = int(r.range(2, 10));
      for (int i = 0; i < k; ++i) {
        const i64 x = r.range(G / 8, 7 * G / 8), y = r.range(G / 8, 7 * G / 8);
        P.ring(star_ring(r, int(r.range(3, max_verts / 2 + 3)), x, y, 1, G / 8, r.chance(80) ? false : true), -1);
      }
      break;
    }
    case 3: {  // random mess: random rings, random roles
      const int k = int(r.range(1, 10));
      for (int i = 0; i < k; ++i)
        P.ring(random_ring(r, int(r.range(1, max_verts / 8 + 3)), 0, 0, G, G), i == 0 || r.chance(40) ? -1 : 0);
      break;
    }
    case 4: {  // touching strips: holes sharing sides (collinear overlaps, not apart)
      const int k = int(r.range(2, 9)), sub = int(r.range(1, max_verts / 4 > 1 ? max_verts / 4 : 1));
      P.ring(rect_ring(G / 16, G / 16, 15 * G / 16, 15 * G / 16, sub, true), -1);
      const i64 w = (12 * G / 16) / k;
      for (int i = 0; i < k; ++i)
        P.ring(rect_ring(2 * G / 16 + i * w, 2 * G / 16, 2 * G / 16 + (i + 1) * w, 14 * G / 16, sub, false), 0);
      break;
    }
    case 5: {  // degenerate rings: repeats, closing repeat, spikes, collinear, fewer than 3
      std::vector<Pt> a = star_ring(r, int(r.range(3, 12)), c, c, G / 8, G / 3, false);
      if (r.chance(50)) a.push_back(a[0]);                   // closing repeat
      a.insert(a.begin() + 1, a[1]);                         // consecutive repeat
      P.ring(a, -1);
      P.ring({{10, 10}, {20, 10}, {10, 10}}, -1);             // fewer than 3 distinct
      P.ring({{5, 5}, {10, 10}, {15, 15}, {20, 20}}, -1);     // zero area
      P.ring({{100, 100}, {200, 100}, {150, 100}, {150, 200}}, -1);  // spike (backtrack along y = 100)
      break;
    }
    case 6: {  // coordinates out of range and huge (int128 validation path when |coord| >= 2^28)
      const i64 big = (i64(1) << 29) + r.range(0, 1000);
      P.ring({{-big, -5}, {big, 7}, {big, big}, {-3, big}}, -1);
      P.ring(star_ring(r, int(r.range(3, 20)), 0, 0, 10, 1 << 20, true), 0);
      if (r.chance(50)) {  // near the int32 limits: products of differences overflow int64
        const i64 M = 2147483647 - r.range(0, 5), m = -M - 1 + r.range(0, 5);
        P.ring({{m, m}, {M, m + 3}, {M - 7, M}, {m + 2, M - 1}, {0, 0}, {M, M - 2}}, -1);
        P.ring({{m + 1, 0}, {M, 1}, {M - 1, 2}, {m + 1, 1}}, -1);
      }
      break;
    }
    case 7: {  // pinches: two rings sharing one vertex, a ring touching itself
      P.ring({{100, 100}, {300, 100}, {300, 300}, {100, 300}}, -1);
      P.ring({{300, 300}, {500, 300}, {500, 500}, {300, 500}}, -1);
      P.ring({{600, 600}, {800, 600}, {700, 700}, {800, 800}, {600, 800}, {700, 700}}, -1);
      break;
    }
    case 8: {  // edges on grid seams: many equal x / y values
      const int k = int(r.range(1, 6));
      for (int i = 0; i < k; ++i) {
        const i64 x = r.range(0, 8) * G / 8, y = r.range(0, 8) * G / 8;
        P.ring(rect_ring(x, y, x + G / 8, y + G / 8, int(r.range(1, 4)), true), -1);
      }
      break;
    }
    default: {  // a big star with nested sibling holes
      P.ring(star_ring(r, nv(), c, c, G / 3, G / 2 - 1, false), -1);
      P.ring(star_ring(r, int(r.range(3, max_verts / 4 + 3)), c, c, G / 20, G / 6, true), 0);
      P.ring(star_ring(r, int(r.range(3, max_verts / 4 + 3)), c, c, 1, G / 25, true), 0);
      break;
    }
  }
  P.end_shape();
}

// The ShapeTotals worst case: an outer square and k holes as vertical strips sharing their long sides (every
// neighbour pair touches, so the witness search scans all 2n witnesses against n edges), n vertices per ring.
inline void strips_shape(Polys& P, i64 G, int k, int n) {
  const int sub = n / 4 > 0 ? n / 4 : 1;
  P.ring(rect_ring(G / 16, G / 16, 15 * G / 16, 15 * G / 16, sub, true), -1);
  const i64 w = (12 * G / 16) / k;
  for (int i = 0; i < k; ++i)
    P.ring(rect_ring(2 * G / 16 + i * w, 2 * G / 16, 2 * G / 16 + (i + 1) * w, 14 * G / 16, sub, false), 0);
  P.end_shape();
}

// ---------------------------------------------------------------- masks
using k2ref::Mask;

inline Mask blank(i32 H, i32 W) {
  Mask m;
  m.H = H, m.W = W, m.px.assign(size_t(H) * W, 0);
  return m;
}

inline Mask random_mask(Rng& r, i32 H, i32 W) {
  Mask m = blank(H, W);
  const int fam = int(r.range(0, 9));
  auto set = [&](i64 x, i64 y) { if (x >= 0 && y >= 0 && x < W && y < H) m.px[size_t(y) * W + size_t(x)] = 1; };
  switch (fam) {
    case 0: break;                                                             // empty
    case 1: std::fill(m.px.begin(), m.px.end(), 1); break;                     // full
    case 2: set(r.range(0, W - 1), r.range(0, H - 1)); break;                  // one pixel
    case 3: for (i32 y = 0; y < H; ++y) for (i32 x = 0; x < W; ++x) m.px[size_t(y) * W + x] = (x + y) & 1; break;
    case 4: {                                                                  // noise
      const int d = int(r.range(1, 99));
      for (auto& v : m.px) v = r.chance(d);
      break;
    }
    case 5: {                                                                  // lines
      const i64 x = r.range(0, W - 1), y = r.range(0, H - 1);
      for (i32 i = 0; i < H; ++i) set(x, i);
      for (i32 i = 0; i < W; ++i) set(i, y);
      for (i32 i = 0; i < W + H; ++i) set(i, i);
      break;
    }
    default: {                                                                 // union of discs (integer test)
      const int k = int(r.range(1, 6));
      for (int t = 0; t < k; ++t) {
        const i64 cx = r.range(0, W - 1), cy = r.range(0, H - 1), rad = r.range(0, (W < H ? W : H) / 2 + 1);
        for (i32 y = 0; y < H; ++y)
          for (i32 x = 0; x < W; ++x)
            if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= rad * rad) m.px[size_t(y) * W + x] = 1;
      }
      if (r.chance(30)) {  // a hole
        const i64 cx = r.range(0, W - 1), cy = r.range(0, H - 1), rad = r.range(0, 4);
        for (i32 y = 0; y < H; ++y)
          for (i32 x = 0; x < W; ++x)
            if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= rad * rad) m.px[size_t(y) * W + x] = 0;
      }
      break;
    }
  }
  return m;
}

// A mask near `base`: random pixel flips along its border region (a model sample against a rater, say).
inline Mask jitter(Rng& r, const Mask& base, int pct) {
  Mask m = base;
  for (i32 y = 0; y < m.H; ++y)
    for (i32 x = 0; x < m.W; ++x) {
      const bool edge = base.at(x, y) != base.at(x + 1, y) || base.at(x, y) != base.at(x, y + 1);
      if (edge && r.chance(pct)) m.px[size_t(y) * m.W + x] ^= 1;
    }
  return m;
}

struct Packed {  // owns a packed batch
  std::vector<u32> words;
  std::vector<i64> word_start;
  std::vector<i32> H, W;
  explicit Packed(const std::vector<Mask>& ms) { k2ref::pack(ms, words, word_start, H, W); }
  k2::Bits bits() const {
    k2::Bits b;
    b.n = i64(H.size()), b.words = words.data(), b.word_start = word_start.data(), b.H = H.data(), b.W = W.data();
    return b;
  }
};

}  // namespace k2gen
