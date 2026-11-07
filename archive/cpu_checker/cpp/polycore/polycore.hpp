// polycore.hpp -- exact integer polygon core: rasterize, graded validity, signed area.
//
// Header-only C++17, no dependencies. Every predicate is exact: int64 on the normal
// input range, __int128 for huge coordinates (chosen per shape). Nothing is repaired;
// every function only measures. Shared with the CUDA engine later (cuda/).
//
// Spec (the slow exact oracles in tests/oracles/ implement the same rules literally):
//  * A shape is a list of rings of integer vertices on a G x G grid (G = 1000) covering
//    a crop of W x H pixels. Grid point (x, y) maps to pixel coordinates (x W / G, y H / G).
//    Rings are implicitly closed.
//  * Pixel (c, r) is inside iff its centre (c + 1/2, r + 1/2) is inside.
//  * Half-open crossing rule: a ray from the centre towards +x; edge (x0,y0)->(x1,y1) counts
//    iff min(y0,y1) <= py < max(y0,y1) and the crossing x is strictly greater than px.
//    Contribution +1 if the edge goes down the screen (y increasing), -1 if up.
//    Even-odd: winding odd. Nonzero: winding != 0. Centres on left/top edges are inside,
//    on right/bottom edges outside: the GPU "top-left" rule [3].
//  * Outer rings are clockwise on screen (y down): positive shoelace sum
//    sum x_i y_{i+1} - x_{i+1} y_i. Holes are counter-clockwise.
//
// Input caps (check_* throw std::invalid_argument beyond them): at most kMaxRings rings,
// kMaxRingVertices vertices per ring, 1 <= W, H <= kMaxDim, 1 <= G <= kMaxGrid.
//
// References:
//  [1] Scanline rendering with per-row edge buckets and x-sorted crossings:
//      https://en.wikipedia.org/wiki/Scanline_rendering ,
//      https://graphics.stanford.edu/courses/cs248-01/scan/scan1.html
//  [2] Exact integer orientation predicates (no floating point, so no adaptive arithmetic
//      needed): https://www.cs.cmu.edu/~quake/robust.html
//  [3] Top-left fill convention: https://learn.microsoft.com/en-us/windows/win32/direct3d11/d3d10-graphics-programming-guide-rasterizer-stage-rules
//  [4] Sort-and-sweep pair pruning (sweep and prune); Bentley-Ottmann was not used because
//      every intersecting pair must be counted anyway:
//      https://en.wikipedia.org/wiki/Sweep_and_prune ,
//      https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <vector>

#if !defined(__SIZEOF_INT128__)
#error "polycore requires __int128 (gcc, clang or nvcc host)"
#endif

namespace polycore {

using i32 = std::int32_t;
using i64 = std::int64_t;
using u8 = std::uint8_t;
using u64 = std::uint64_t;
using i128 = __int128;

constexpr i32 kMaxDim = 4096;
constexpr i32 kMaxRings = 10;
constexpr i32 kMaxRingVertices = 1000;
constexpr i64 kMaxGrid = i64(1) << 30;

// ---------------------------------------------------------------- input view
// Vertices are flat int32 pairs x0,y0,x1,y1,...; ring r owns vertices [off[r], off[r+1]).
struct Rings {
  const i32* xy = nullptr;
  const i32* off = nullptr;
  i32 nrings = 0;
  i32 begin(i32 r) const { return off[r]; }
  i32 size(i32 r) const { return off[r + 1] - off[r]; }
  i32 x(i32 v) const { return xy[2 * v]; }
  i32 y(i32 v) const { return xy[2 * v + 1]; }
};

// n_vertices = number of (x, y) pairs available in rs.xy.
inline void check_rings(const Rings& rs, i64 n_vertices) {
  if (rs.nrings < 0 || rs.nrings > kMaxRings) throw std::invalid_argument("at most 10 rings per shape");
  if (rs.nrings == 0) return;
  if (rs.off[0] < 0) throw std::invalid_argument("ring offsets must start at >= 0");
  for (i32 r = 0; r < rs.nrings; ++r) {
    if (rs.off[r + 1] < rs.off[r]) throw std::invalid_argument("ring offsets must be non-decreasing");
    if (rs.size(r) > kMaxRingVertices) throw std::invalid_argument("at most 1000 vertices per ring");
  }
  if (rs.off[rs.nrings] > n_vertices) throw std::invalid_argument("ring offsets exceed the number of vertices");
}

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

namespace detail {

template <class T>
inline T floordiv(T a, T b) {  // b != 0
  T q = a / b;
  if ((a % b != 0) && ((a < 0) != (b < 0))) --q;
  return q;
}

inline i64 iabs(i64 v) { return v < 0 ? -v : v; }

inline i64 max_abs_coord(const Rings& rs) {
  i64 m = 0;
  if (rs.nrings <= 0) return 0;
  for (i32 v = rs.off[0]; v < rs.off[rs.nrings]; ++v) m = std::max(m, std::max(iabs(rs.x(v)), iabs(rs.y(v))));
  return m;
}

inline i64 saturate(i128 v) {
  const i128 hi = std::numeric_limits<i64>::max(), lo = std::numeric_limits<i64>::min();
  return v > hi ? i64(hi) : (v < lo ? i64(lo) : i64(v));
}

}  // namespace detail

// Twice the signed shoelace area of ring r (exact). > 0 means clockwise on screen.
inline i128 signed_area2(const Rings& rs, i32 r) {
  const i32 b = rs.begin(r), n = rs.size(r);
  i128 acc = 0;
  for (i32 k = 0; k < n; ++k) {
    const i32 v0 = b + k, v1 = b + (k + 1 == n ? 0 : k + 1);
    acc += i128(rs.x(v0)) * rs.y(v1) - i128(rs.x(v1)) * rs.y(v0);
  }
  return acc;
}

// ================================================================= rasterize
// Each edge crossing a pixel-centre row becomes one key (row, col, ring, sign), where col is
// the number of centres of that row strictly left of the crossing (clamped to [0, W]): the
// edge adds its sign to the winding of exactly those centres. Keys are sorted, and each row
// is walked left to right as spans of constant winding [1]. Memory is O(crossings), never
// O(W H); rows are processed in bands so that a pathological shape stays bounded.
struct RasterScratch {
  std::vector<u64> keys;
};

struct RasterOut {           // optional outputs, any pointer may be null
  u8* evenodd = nullptr;     // H*W, 0/1
  u8* nonzero = nullptr;     // H*W, 0/1
  i32* ring_wmin = nullptr;  // per ring: min(0, min of the ring's own winding over centres)
  i32* ring_wmax = nullptr;  // per ring: max(0, max ...)
};

struct RasterStats {
  i64 area_evenodd = 0;      // pixels with odd winding
  i64 area_nonzero = 0;      // pixels with nonzero winding
  i64 px_multi = 0;          // pixels with |winding| >= 2 (self-overlap of one orientation)
  i32 rings_mixed_sign = 0;  // rings whose own winding takes both signs (bow-tie, figure-8)
  i64 gt_area = 0, inter_evenodd = 0, inter_nonzero = 0;  // filled only with a gt mask
};

namespace detail {

constexpr i32 kBandRows = 512;  // bounds keys per band to 10 * 1000 * 512
constexpr int kColShift = 5, kRowShift = 18;  // key = row << 18 | col << 5 | ring << 1 | (sign > 0)

// Keys of all edge crossings with centre rows in [b0, b1).
template <class T>
inline void collect_keys(const Rings& rs, i32 W, i32 H, i64 G, i32 b0, i32 b1, std::vector<u64>& keys) {
  // Scaled frame: x * 2W, y * 2H; centre (c, r) is ((2c+1) G, (2r+1) G). All integers.
  const T sx = T(2) * W, sy = T(2) * H, g = T(G), g2 = T(2) * G;
  for (i32 ring = 0; ring < rs.nrings; ++ring) {
    const i32 b = rs.begin(ring), n = rs.size(ring);
    for (i32 k = 0; k < n; ++k) {
      const i32 v0 = b + k, v1 = b + (k + 1 == n ? 0 : k + 1);
      const T ay = T(rs.y(v0)) * sy, by = T(rs.y(v1)) * sy;
      if (ay == by) continue;  // horizontal edges never count
      const T ax = T(rs.x(v0)) * sx, bx = T(rs.x(v1)) * sx;
      const u64 up = by > ay;  // edge goes down the screen: +1
      const T x0 = ay < by ? ax : bx, y0 = ay < by ? ay : by;
      const T dx = (ay < by ? bx : ax) - x0, dy = (ay < by ? by : ay) - y0;  // dy > 0
      // rows with y0 <= (2r+1) G < y0 + dy (half-open)
      T rlo = -floordiv<T>(g - y0, g2), rhi = floordiv<T>(y0 + dy - 1 - g, g2);
      if (rlo < T(b0)) rlo = T(b0);
      if (rhi > T(b1 - 1)) rhi = T(b1 - 1);
      if (rlo > rhi) continue;
      const T gdy = g * dy;
      for (i32 r = i32(rlo); r <= i32(rhi); ++r) {
        const T num = x0 * dy + (T(2 * r + 1) * g - y0) * dx;  // crossing_x * dy
        // centre c counts iff (2c+1) G dy < num; q = max k with k G dy < num
        const T q = floordiv<T>(num - 1, gdy);
        const T cnt = floordiv<T>(q - 1, T(2)) + 1;
        const u64 col = cnt < 0 ? 0 : (cnt > T(W) ? u64(W) : u64(cnt));
        keys.push_back(u64(r) << kRowShift | col << kColShift | u64(ring) << 1 | up);
      }
    }
  }
}

// Walk sorted keys row by row as spans [a, b) of constant winding.
inline void walk_spans(const std::vector<u64>& keys, i32 W, i32 R, const RasterOut& out, const u8* gt,
                       RasterStats& st, i32* wmin, i32* wmax) {
  const size_t n = keys.size();
  size_t i = 0;
  while (i < n) {
    const u64 row = keys[i] >> kRowShift;
    const size_t base = size_t(row) * size_t(W);
    i32 tot = 0, wr[kMaxRings] = {};
    while (i < n && (keys[i] >> kRowShift) == row) {
      const u64 pos = keys[i] >> kColShift;  // row and col
      const i32 a = i32(pos & 0x1FFF);
      // winding(c) = -(sum of signs with col <= c): a closed ring's signs sum to 0 per row
      for (; i < n && (keys[i] >> kColShift) == pos; ++i) {
        const i32 s = (keys[i] & 1) ? 1 : -1;
        tot -= s;
        wr[(keys[i] >> 1) & 0xF] -= s;
      }
      const i32 b = (i < n && (keys[i] >> kRowShift) == row) ? i32((keys[i] >> kColShift) & 0x1FFF) : W;
      if (b <= a) continue;
      for (i32 r = 0; r < R; ++r) {
        wmin[r] = std::min(wmin[r], wr[r]);
        wmax[r] = std::max(wmax[r], wr[r]);
      }
      if (tot == 0) continue;
      const i64 len = b - a;
      const bool odd = tot & 1;
      st.area_nonzero += len;
      st.area_evenodd += odd ? len : 0;
      st.px_multi += (tot >= 2 || tot <= -2) ? len : 0;
      if (out.nonzero) std::memset(out.nonzero + base + a, 1, size_t(len));
      if (out.evenodd && odd) std::memset(out.evenodd + base + a, 1, size_t(len));
      if (gt) {
        i64 cnt = 0;
        for (i32 c = a; c < b; ++c) cnt += gt[base + c] != 0;
        st.inter_nonzero += cnt;
        st.inter_evenodd += odd ? cnt : 0;
      }
    }
  }
}

}  // namespace detail

// True iff int64 arithmetic is exact for this shape on this crop (all scaled coordinates,
// centres and the products formed in collect_keys stay below 2^62).
inline bool raster_fits_i64(const Rings& rs, i32 W, i32 H, i64 G) {
  const i64 dim = std::max(W, H), lim = i64(1) << 30;
  return detail::max_abs_coord(rs) * 2 * dim < lim && 2 * dim * G < lim;
}

// Rasterize all rings together (fill rules apply to the union of rings). gt (optional) is
// an H*W mask, used only to count intersections. Inputs must pass check_rings/check_dims/
// check_grid.
inline RasterStats rasterize(const Rings& rs, i32 W, i32 H, i64 G, RasterScratch& sc, const RasterOut& out,
                             const u8* gt = nullptr) {
  RasterStats st;
  const size_t npx = size_t(W) * size_t(H);
  if (out.evenodd) std::memset(out.evenodd, 0, npx);
  if (out.nonzero) std::memset(out.nonzero, 0, npx);
  i32 wmin[kMaxRings] = {}, wmax[kMaxRings] = {};
  const bool small = raster_fits_i64(rs, W, H, G);
  for (i32 b0 = 0; b0 < H; b0 += detail::kBandRows) {
    const i32 b1 = std::min(H, b0 + detail::kBandRows);
    sc.keys.clear();
    if (small) detail::collect_keys<i64>(rs, W, H, G, b0, b1, sc.keys);
    else detail::collect_keys<i128>(rs, W, H, G, b0, b1, sc.keys);
    std::sort(sc.keys.begin(), sc.keys.end());
    detail::walk_spans(sc.keys, W, rs.nrings, out, gt, st, wmin, wmax);
  }
  for (i32 r = 0; r < rs.nrings; ++r) {
    if (out.ring_wmin) out.ring_wmin[r] = wmin[r];
    if (out.ring_wmax) out.ring_wmax[r] = wmax[r];
    st.rings_mixed_sign += wmin[r] < 0 && wmax[r] > 0;
  }
  if (gt) {
    i64 a = 0;
    for (size_t p = 0; p < npx; ++p) a += gt[p] != 0;
    st.gt_area = a;
  }
  return st;
}
// ================================================================== validate
enum Defect : int {
  kFewerThan3Points = 0,  // < 3 distinct vertices after dropping repeats (no other self label then)
  kZeroArea,              // shoelace sum == 0
  kSelfCrossing,          // # non-adjacent edge pairs of the ring that cross properly
  kVertexPinch,           // # non-adjacent pairs touching at a point (vertex=vertex or vertex on edge)
  kCollinearOverlap,      // # non-adjacent pairs overlapping along a segment of positive length
  kSpike,                 // # adjacent pairs folding back on each other (collinear backtrack)
  kRingRingCrossing,      // # edge pairs (this ring x other rings) crossing properly
  kRingRingTouching,      // # edge pairs (this ring x other rings) touching / overlapping
  kHoleOutsideOuter,      // hole with a witness point strictly outside its outer ring (a hole that
                          // crosses its outer gets this flag and ring_ring_crossing)
  kNestedHoles,           // # sibling holes this hole overlaps / is nested with
  kNestedOuters,          // # other outer parts this outer overlaps / is nested with
  kWrongWinding,          // outer with net area < 0 or hole with net area > 0 (net shoelace area
                          // of the cleaned ring, so a bow-tie is judged by its larger lobe)
  kOutOfRange,            // # vertices with x or y outside [0, G]
  kNumDefects
};

inline constexpr const char* kDefectNames[kNumDefects] = {
    "fewer_than_3_points", "zero_area",          "self_crossing",
    "vertex_pinch",        "collinear_overlap",  "spike",
    "ring_ring_crossing",  "ring_ring_touching", "hole_outside_outer",
    "nested_holes",        "nested_outers",      "wrong_winding",
    "out_of_range"};

struct ValidateResult {
  i32 nrings = 0;
  std::vector<i32> counts;  // nrings x kNumDefects
  std::array<i64, kNumDefects> total{};  // ring-ring pair classes counted once per pair
  std::vector<i64> area2;   // signed twice-area (saturated to int64)
  std::vector<i32> n_clean; // vertices after dropping consecutive / closing repeats
  std::vector<i32> pair_cross, pair_touch;  // nrings x nrings, symmetric, diag 0
  void reset(i32 R) {
    nrings = R;
    counts.assign(size_t(R) * kNumDefects, 0);
    total.fill(0);
    area2.assign(R, 0);
    n_clean.assign(R, 0);
    pair_cross.assign(size_t(R) * R, 0);
    pair_touch.assign(size_t(R) * R, 0);
  }
  i32& at(i32 r, Defect d) { return counts[size_t(r) * kNumDefects + d]; }
};

struct ValidateScratch {
  struct Edge {
    i64 ax, ay, bx, by, xmin, xmax, ymin, ymax;
    i32 ring, idx, n;
  };
  std::vector<i64> px, py;        // clean vertices of all rings
  std::vector<i32> start, len;    // clean ring r = [start[r], start[r]+len[r])
  std::vector<u8> fewer;
  std::vector<i64> bx0, bx1, by0, by1;  // clean ring bounding boxes
  std::vector<Edge> edges;
  std::vector<i32> order;
  std::vector<i32> children;      // hole lists per outer (flattened)
  std::vector<i32> child_start;
};

namespace detail {

template <class T>
inline int orient(i64 ax, i64 ay, i64 bx, i64 by, i64 cx, i64 cy) {
  const T v = T(bx - ax) * T(cy - ay) - T(by - ay) * T(cx - ax);
  return (v > 0) - (v < 0);
}

inline bool in_box(i64 ax, i64 ay, i64 bx, i64 by, i64 px, i64 py) {
  return std::min(ax, bx) <= px && px <= std::max(ax, bx) && std::min(ay, by) <= py &&
         py <= std::max(ay, by);
}

// p on closed segment ab (exact)
template <class T>
inline bool on_segment(i64 ax, i64 ay, i64 bx, i64 by, i64 px, i64 py) {
  return in_box(ax, ay, bx, by, px, py) && orient<T>(ax, ay, bx, by, px, py) == 0;
}

enum PairKind { kNone = 0, kCross, kPinch, kOverlap, kSpikeKind };

// Classify two edges. Same-ring adjacency follows oracle.ring_problems exactly.
template <class T>
inline PairKind classify(const ValidateScratch::Edge& e, const ValidateScratch::Edge& f) {
  if (e.ring == f.ring) {
    const i32 n = e.n;
    bool e_then_f = (f.idx == (e.idx + 1 == n ? 0 : e.idx + 1));
    bool f_then_e = (e.idx == (f.idx + 1 == n ? 0 : f.idx + 1));
    if (e_then_f || f_then_e) {
      // shared vertex S, other endpoints P (first edge) and Q (second edge)
      const auto& a = e_then_f ? e : f;  // edge ending at S
      const auto& b = e_then_f ? f : e;  // edge starting at S
      const i64 sx = a.bx, sy = a.by, px = a.ax, py = a.ay, qx = b.bx, qy = b.by;
      if (orient<T>(sx, sy, px, py, qx, qy) == 0 &&
          (in_box(sx, sy, px, py, qx, qy) || in_box(sx, sy, qx, qy, px, py)))
        return kSpikeKind;
      return kNone;
    }
  }
  const int o1 = orient<T>(e.ax, e.ay, e.bx, e.by, f.ax, f.ay);
  const int o2 = orient<T>(e.ax, e.ay, e.bx, e.by, f.bx, f.by);
  const int o3 = orient<T>(f.ax, f.ay, f.bx, f.by, e.ax, e.ay);
  const int o4 = orient<T>(f.ax, f.ay, f.bx, f.by, e.bx, e.by);
  if (o1 * o2 < 0 && o3 * o4 < 0) return kCross;
  const bool touch = (o1 == 0 && in_box(e.ax, e.ay, e.bx, e.by, f.ax, f.ay)) ||
                     (o2 == 0 && in_box(e.ax, e.ay, e.bx, e.by, f.bx, f.by)) ||
                     (o3 == 0 && in_box(f.ax, f.ay, f.bx, f.by, e.ax, e.ay)) ||
                     (o4 == 0 && in_box(f.ax, f.ay, f.bx, f.by, e.bx, e.by));
  if (!touch) return kNone;
  if (o1 == 0 && o2 == 0) {  // collinear: overlap of positive length?
    const bool use_x = e.ax != e.bx;
    const i64 lo = use_x ? std::max(e.xmin, f.xmin) : std::max(e.ymin, f.ymin);
    const i64 hi = use_x ? std::min(e.xmax, f.xmax) : std::min(e.ymax, f.ymax);
    if (hi > lo) return kOverlap;
  }
  return kPinch;
}

inline bool boxes_meet(const ValidateScratch::Edge& e, const ValidateScratch::Edge& f) {
  return e.xmin <= f.xmax && f.xmin <= e.xmax && e.ymin <= f.ymax && f.ymin <= e.ymax;
}

template <class T>
inline void record_pair(const ValidateScratch::Edge& e, const ValidateScratch::Edge& f,
                        ValidateResult& out) {
  const PairKind k = classify<T>(e, f);
  if (k == kNone) return;
  const i32 R = out.nrings;
  if (e.ring == f.ring) {
    Defect d = k == kCross ? kSelfCrossing
             : k == kPinch ? kVertexPinch
             : k == kOverlap ? kCollinearOverlap : kSpike;
    out.at(e.ring, d) += 1;
  } else {
    i32* m = (k == kCross) ? out.pair_cross.data() : out.pair_touch.data();
    m[size_t(e.ring) * R + f.ring] += 1;
    m[size_t(f.ring) * R + e.ring] += 1;
  }
}

// Locate a point given in doubled coordinates against clean ring r:
// -1 strictly outside, 0 on the boundary, +1 strictly inside (even-odd).
template <class T>
inline int locate2(const ValidateScratch& s, i32 r, i64 qx, i64 qy) {
  if (qx < 2 * s.bx0[r] || qx > 2 * s.bx1[r] || qy < 2 * s.by0[r] || qy > 2 * s.by1[r]) return -1;
  const i32 b = s.start[r], n = s.len[r];
  bool in = false;
  for (i32 k = 0; k < n; ++k) {
    const i32 i0 = b + k, i1 = b + (k + 1 == n ? 0 : k + 1);
    const i64 ax = 2 * s.px[i0], ay = 2 * s.py[i0], bx = 2 * s.px[i1], by = 2 * s.py[i1];
    if (on_segment<T>(ax, ay, bx, by, qx, qy)) return 0;
    if ((ay <= qy && qy < by) || (by <= qy && qy < ay)) {
      const T num = T(ax - qx) * T(by - ay) + T(qy - ay) * T(bx - ax);
      if (by > ay ? num > 0 : num < 0) in = !in;
    }
  }
  return in ? 1 : -1;
}

// Region of an outer ring with its holes (even-odd over all of them).
template <class T>
inline int locate2_region(const ValidateScratch& s, i32 outer, i64 qx, i64 qy) {
  int loc = locate2<T>(s, outer, qx, qy);
  if (loc == 0) return 0;
  bool in = loc > 0;
  for (i32 k = s.child_start[outer]; k < s.child_start[outer + 1]; ++k) {
    const i32 h = s.children[k];
    if (s.fewer[h]) continue;
    const int lh = locate2<T>(s, h, qx, qy);
    if (lh == 0) return 0;
    if (lh > 0) in = !in;
  }
  return in ? 1 : -1;
}

// Does some witness point of ring a (vertex or edge midpoint) satisfy loc(w) == want?
// single: ring a neither crosses nor touches the boundary being tested, so its whole
// boundary lies in one component of the complement and one vertex decides exactly.
template <class T, class Loc>
inline bool any_witness(const ValidateScratch& s, i32 a, int want, bool single, Loc loc) {
  const i32 b = s.start[a], n = s.len[a];
  if (single) return n > 0 && loc(2 * s.px[b], 2 * s.py[b]) == want;
  for (i32 k = 0; k < n; ++k) {
    const i32 i0 = b + k, i1 = b + (k + 1 == n ? 0 : k + 1);
    if (loc(2 * s.px[i0], 2 * s.py[i0]) == want) return true;
    if (loc(s.px[i0] + s.px[i1], s.py[i0] + s.py[i1]) == want) return true;
  }
  return false;
}

template <class T>
inline void validate_geometry(const Rings& rs, const i32* roles, ValidateResult& out,
                              ValidateScratch& s) {
  const i32 R = rs.nrings;
  // --- all edges of all usable rings
  s.edges.clear();
  for (i32 r = 0; r < R; ++r) {
    if (s.fewer[r]) continue;
    const i32 b = s.start[r], n = s.len[r];
    for (i32 k = 0; k < n; ++k) {
      const i32 i0 = b + k, i1 = b + (k + 1 == n ? 0 : k + 1);
      ValidateScratch::Edge e;
      e.ax = s.px[i0]; e.ay = s.py[i0]; e.bx = s.px[i1]; e.by = s.py[i1];
      e.xmin = std::min(e.ax, e.bx); e.xmax = std::max(e.ax, e.bx);
      e.ymin = std::min(e.ay, e.by); e.ymax = std::max(e.ay, e.by);
      e.ring = r; e.idx = k; e.n = n;
      s.edges.push_back(e);
    }
  }
  const i32 E = i32(s.edges.size());
  const auto& ed = s.edges;
  // sort-and-sweep on x, test y-overlap [4]
  {
    s.order.resize(E);
    for (i32 i = 0; i < E; ++i) s.order[i] = i;
    std::sort(s.order.begin(), s.order.end(), [&](i32 a, i32 b) {
      return ed[a].xmin != ed[b].xmin ? ed[a].xmin < ed[b].xmin : a < b;
    });
    for (i32 ii = 0; ii < E; ++ii) {
      const auto& e = ed[s.order[ii]];
      for (i32 jj = ii + 1; jj < E; ++jj) {
        const auto& f = ed[s.order[jj]];
        if (f.xmin > e.xmax) break;
        if (f.ymin <= e.ymax && e.ymin <= f.ymax) record_pair<T>(e, f, out);
      }
    }
  }
  // --- ring-ring counts per ring and totals
  for (i32 r = 0; r < R; ++r)
    for (i32 q = 0; q < R; ++q) {
      out.at(r, kRingRingCrossing) += out.pair_cross[size_t(r) * R + q];
      out.at(r, kRingRingTouching) += out.pair_touch[size_t(r) * R + q];
    }
  // --- roles: children lists (holes of each outer)
  s.child_start.assign(R + 1, 0);
  s.children.clear();
  if (roles) {
    for (i32 r = 0; r < R; ++r)
      if (roles[r] >= 0) s.child_start[roles[r] + 1] += 1;
    for (i32 r = 0; r < R; ++r) s.child_start[r + 1] += s.child_start[r];
    s.children.assign(s.child_start[R], 0);
    std::vector<i32> fill(s.child_start.begin(), s.child_start.end() - 1);
    for (i32 r = 0; r < R; ++r)
      if (roles[r] >= 0) s.children[fill[roles[r]]++] = r;
  }
  auto is_hole = [&](i32 r) { return roles && roles[r] >= 0; };
  // --- winding
  for (i32 r = 0; r < R; ++r) {
    if (s.fewer[r]) continue;
    const i64 a = out.area2[r];
    if ((!is_hole(r) && a < 0) || (is_hole(r) && a > 0)) out.at(r, kWrongWinding) = 1;
  }
  // --- containment / nesting via exact witness points
  auto apart = [&](i32 a, i32 b) {
    return out.pair_cross[size_t(a) * R + b] == 0 && out.pair_touch[size_t(a) * R + b] == 0;
  };
  auto apart_region = [&](i32 a, i32 o) {  // a vs outer o and all of o's holes
    if (!apart(a, o)) return false;
    for (i32 k = s.child_start[o]; k < s.child_start[o + 1]; ++k)
      if (!apart(a, s.children[k])) return false;
    return true;
  };
  for (i32 h = 0; h < R; ++h) {
    if (!is_hole(h) || s.fewer[h]) continue;
    const i32 p = roles[h];
    if (s.fewer[p]) continue;
    if (any_witness<T>(s, h, -1, apart(h, p), [&](i64 x, i64 y) { return locate2<T>(s, p, x, y); }))
      out.at(h, kHoleOutsideOuter) = 1;
  }
  for (i32 a = 0; a < R; ++a) {
    if (s.fewer[a]) continue;
    for (i32 b = a + 1; b < R; ++b) {
      if (s.fewer[b]) continue;
      const bool both_holes = is_hole(a) && is_hole(b) && roles[a] == roles[b];
      const bool both_outers = !is_hole(a) && !is_hole(b);
      bool nested = false;
      if (both_holes) {
        const bool single = apart(a, b);
        nested = any_witness<T>(s, b, 1, single, [&](i64 x, i64 y) { return locate2<T>(s, a, x, y); }) ||
                 any_witness<T>(s, a, 1, single, [&](i64 x, i64 y) { return locate2<T>(s, b, x, y); });
        if (nested) { out.at(a, kNestedHoles) += 1; out.at(b, kNestedHoles) += 1; out.total[kNestedHoles] += 1; }
      } else if (both_outers) {
        nested = any_witness<T>(s, b, 1, apart_region(b, a), [&](i64 x, i64 y) { return locate2_region<T>(s, a, x, y); }) ||
                 any_witness<T>(s, a, 1, apart_region(a, b), [&](i64 x, i64 y) { return locate2_region<T>(s, b, x, y); });
        if (nested) { out.at(a, kNestedOuters) += 1; out.at(b, kNestedOuters) += 1; out.total[kNestedOuters] += 1; }
      }
    }
  }
}

}  // namespace detail

// Graded legality of a shape. roles: nullptr => every ring is an outer ring;
// otherwise roles[r] == -1 for an outer ring and roles[r] == p (an outer ring
// index) for a hole of ring p (see check_roles). G is the grid size
// used for out_of_range. Never repairs anything.
inline void validate(const Rings& rs, const i32* roles, i64 G, ValidateResult& out,
                     ValidateScratch& s) {
  const i32 R = rs.nrings;
  out.reset(R);
  s.px.clear(); s.py.clear();
  s.start.assign(R, 0); s.len.assign(R, 0); s.fewer.assign(R, 0);
  s.bx0.assign(R, 0); s.bx1.assign(R, 0); s.by0.assign(R, 0); s.by1.assign(R, 0);
  for (i32 r = 0; r < R; ++r) {
    const i32 b = rs.begin(r), n = rs.size(r);
    const i32 st = i32(s.px.size());
    s.start[r] = st;
    i32 oor = 0;
    for (i32 k = 0; k < n; ++k) {
      const i64 x = rs.x(b + k), y = rs.y(b + k);
      oor += (x < 0 || x > G || y < 0 || y > G);
      if (k == 0 || x != rs.x(b + k - 1) || y != rs.y(b + k - 1)) { s.px.push_back(x); s.py.push_back(y); }
    }
    i32 m = i32(s.px.size()) - st;
    if (m > 1 && s.px[st] == s.px[st + m - 1] && s.py[st] == s.py[st + m - 1]) {
      s.px.pop_back(); s.py.pop_back(); --m;
    }
    s.len[r] = m;
    out.n_clean[r] = m;
    if (m > 0) {
      s.bx0[r] = *std::min_element(s.px.begin() + st, s.px.end());
      s.bx1[r] = *std::max_element(s.px.begin() + st, s.px.end());
      s.by0[r] = *std::min_element(s.py.begin() + st, s.py.end());
      s.by1[r] = *std::max_element(s.py.begin() + st, s.py.end());
    }
    out.at(r, kOutOfRange) = oor;
    // fewer than 3 distinct points?
    i32 distinct = m > 0 ? 1 : 0, second = -1;
    for (i32 k = 1; k < m && distinct < 3; ++k) {
      const i64 x = s.px[st + k], y = s.py[st + k];
      if (x == s.px[st] && y == s.py[st]) continue;
      if (second < 0) { second = st + k; distinct = 2; continue; }
      if (x == s.px[second] && y == s.py[second]) continue;
      distinct = 3;
    }
    i128 a2 = 0;
    for (i32 k = 0; k < m; ++k) {
      const i32 i0 = st + k, i1 = st + (k + 1 == m ? 0 : k + 1);
      a2 += i128(s.px[i0]) * s.py[i1] - i128(s.px[i1]) * s.py[i0];
    }
    out.area2[r] = detail::saturate(a2);
    if (distinct < 3) {
      s.fewer[r] = 1;
      out.at(r, kFewerThan3Points) = 1;
    } else if (a2 == 0) {
      out.at(r, kZeroArea) = 1;
    }
  }
  if (detail::max_abs_coord(rs) < (i64(1) << 28)) detail::validate_geometry<i64>(rs, roles, out, s);
  else detail::validate_geometry<i128>(rs, roles, out, s);
  for (i32 r = 0; r < R; ++r)
    for (int d = 0; d < kNumDefects; ++d)
      if (d != kRingRingCrossing && d != kRingRingTouching && d != kNestedHoles && d != kNestedOuters)
        out.total[d] += out.at(r, Defect(d));
  for (i32 r = 0; r < R; ++r)
    for (i32 q = r + 1; q < R; ++q) {
      out.total[kRingRingCrossing] += out.pair_cross[size_t(r) * R + q];
      out.total[kRingRingTouching] += out.pair_touch[size_t(r) * R + q];
    }
}

// ===================================================================== batch
struct Rollout {
  Rings rings;
  const i32* roles = nullptr;
  i32 W = 0, H = 0;
  const u8* gt = nullptr;  // optional H*W mask
};

struct RolloutScore {
  RasterStats raster;
  std::array<i64, kNumDefects> defects{};
};

struct Workspace {
  RasterScratch raster;
  ValidateScratch val;
  ValidateResult res;
};

inline void score_one(const Rollout& ro, i64 G, Workspace& ws, RolloutScore& out) {
  out.raster = rasterize(ro.rings, ro.W, ro.H, G, ws.raster, RasterOut{}, ro.gt);
  validate(ro.rings, ro.roles, G, ws.res, ws.val);
  out.defects = ws.res.total;
}

// Score n rollouts on `threads` worker threads (<= 1: run inline). Results are
// written by index, so they do not depend on the thread count or schedule.
// An exception escaping a std::thread calls std::terminate [5], so each worker
// catches it into its own std::exception_ptr and stops the others; the first one
// (by worker index) is rethrown on the calling thread after every worker has
// joined (join synchronizes with the worker's completion [6]). If a thread cannot
// be started (std::system_error), the threads already running finish the batch.
//  [5] https://en.cppreference.com/w/cpp/thread/thread/thread
//  [6] https://en.cppreference.com/w/cpp/error/exception_ptr
inline void score_batch(const Rollout* in, i64 n, i64 G, int threads, RolloutScore* out) {
  if (threads <= 1 || n <= 1) {
    Workspace ws;
    for (i64 i = 0; i < n; ++i) score_one(in[i], G, ws, out[i]);
    return;
  }
  const int nt = int(std::min<i64>(threads, n));
  std::atomic<i64> next{0};
  std::vector<std::exception_ptr> err(static_cast<size_t>(nt));
  auto work = [&](int t) {
    try {
      Workspace ws;
      for (i64 i; (i = next.fetch_add(1, std::memory_order_relaxed)) < n;) score_one(in[i], G, ws, out[i]);
    } catch (...) {
      err[size_t(t)] = std::current_exception();
      next.store(n, std::memory_order_relaxed);  // the other workers stop at their next rollout
    }
  };
  std::vector<std::thread> pool;
  pool.reserve(size_t(nt - 1));
  try {
    for (int t = 1; t < nt; ++t) pool.emplace_back(work, t);
  } catch (const std::system_error&) {
  }
  work(0);
  for (auto& th : pool) th.join();
  for (const auto& e : err)
    if (e) std::rethrow_exception(e);
}

}  // namespace polycore
