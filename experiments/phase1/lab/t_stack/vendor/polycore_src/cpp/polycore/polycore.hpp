// polycore.hpp -- exact integer polygon core for the geometry reward checker.
//
// Header-only, C++17, no dependencies. All predicates are exact: int64 on the
// normal input range, __int128 for garbage-sized coordinates (dispatched per
// shape). Nothing is ever repaired; every function only measures.
//
// Conventions (identical to phase1_lab/common/oracle.py, which is the spec):
//  * A shape is a list of rings of integer vertices on a G x G grid (G = 1000)
//    covering a crop of W x H pixels. Grid point (x, y) maps to pixel
//    coordinates (x * W / G, y * H / G). Rings are implicitly closed.
//  * Pixel (c, r) is inside iff its centre (c + 1/2, r + 1/2) is inside.
//  * Half-open crossing rule: ray from the centre towards +x; edge
//    (x0,y0)->(x1,y1) counts iff min(y0,y1) <= py < max(y0,y1) and the crossing
//    x is strictly greater than px. Winding contribution +1 if the edge goes
//    down the screen (y increasing), -1 if up. Even-odd: winding odd. Nonzero:
//    winding != 0. ("Top-left" tie rule: centres on left/top edges are in.)
//  * Orientation: outer rings are clockwise on screen (y down), i.e. positive
//    shoelace sum  sum x_i*y_{i+1} - x_{i+1}*y_i  ; holes counter-clockwise.
//
// Exact frame used by the rasterizer: scale x by 2W and y by 2H, then pixel
// centres are ((2c+1) G, (2r+1) G) -- every quantity is an integer.
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <limits>
#include <thread>
#include <vector>

#if !defined(__SIZEOF_INT128__)
#error "polycore requires __int128 (gcc, clang or nvcc host)"
#endif

namespace polycore {

using i32 = std::int32_t;
using i64 = std::int64_t;
using u8 = std::uint8_t;
using i128 = __int128;

constexpr i32 kMaxDim = 1 << 12;          // max W, H (scratch is O(W*H) int32)
constexpr i64 kMaxGrid = i64(1) << 30;    // max G
constexpr i32 kMaxRingsValidate = 1024;   // validate() keeps R x R pair matrices

// ---------------------------------------------------------------- input view
// Vertices are flat int32 pairs x0,y0,x1,y1,...; ring r owns vertices
// [off[r], off[r+1]). The caller guarantees off is non-decreasing, off[0] >= 0
// and 2*off[nrings] <= length of xy (the Python binding checks this).
struct Rings {
  const i32* xy = nullptr;
  const i32* off = nullptr;
  i32 nrings = 0;
  i32 begin(i32 r) const { return off[r]; }
  i32 size(i32 r) const { return off[r + 1] - off[r]; }
  i32 x(i32 v) const { return xy[2 * v]; }
  i32 y(i32 v) const { return xy[2 * v + 1]; }
};

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
  const i32 v0 = rs.off[0], v1 = rs.off[rs.nrings];
  for (i32 v = v0; v < v1; ++v) m = std::max(m, std::max(iabs(rs.x(v)), iabs(rs.y(v))));
  return m;
}

inline i64 saturate(i128 v) {
  const i128 hi = std::numeric_limits<i64>::max(), lo = std::numeric_limits<i64>::min();
  return v > hi ? i64(hi) : (v < lo ? i64(lo) : i64(v));
}

}  // namespace detail

// ---------------------------------------------------------------- signed area
// Twice the signed shoelace area (exact). > 0 means clockwise on screen.
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
struct RasterScratch {
  std::vector<i32> diff;  // H x (W+1) signed crossing differences; all-zero between calls
  std::vector<i32> wind;  // H x W   accumulated winding;          all-zero between calls
  void ensure(i32 W, i32 H) {
    const size_t nd = size_t(H) * size_t(W + 1), nw = size_t(H) * size_t(W);
    if (diff.size() < nd) diff.assign(nd, 0);
    if (wind.size() < nw) wind.assign(nw, 0);
  }
};

// Optional outputs: any pointer may be null.
struct RasterOut {
  i32* winding = nullptr;    // H*W
  u8* evenodd = nullptr;     // H*W, 0/1
  u8* nonzero = nullptr;     // H*W, 0/1
  i32* ring_wmin = nullptr;  // nrings: min of the ring's own winding over pixel centres
  i32* ring_wmax = nullptr;  // nrings: max ...
};

struct RasterStats {
  i64 area_evenodd = 0;  // pixels with odd winding
  i64 area_nonzero = 0;  // pixels with nonzero winding
  i64 px_multi = 0;      // pixels with |winding| >= 2 (self-overlap of same orientation)
  i32 rings_mixed_sign = 0;  // rings whose own winding takes both signs (bow-tie, figure-8)
  // filled only when a ground-truth mask is given
  i64 gt_area = 0, inter_evenodd = 0, inter_nonzero = 0;
};

namespace detail {

// Accumulate the crossing differences of one ring into diff (row stride W+1).
// Entry -s at column (number of pixel centres strictly left of the crossing),
// clamped to [0, W]; the matching +s at column 0 is omitted because the signed
// crossings of a closed ring sum to zero on every row.
template <class T>
inline void scan_ring(const Rings& rs, i32 ring, i32 W, i32 H, i64 G, i32* diff, i32& rmin,
                      i32& rmax, i32& cmin, i32& cmax) {
  const i32 b = rs.begin(ring), n = rs.size(ring);
  const T sx = T(2) * W, sy = T(2) * H, g = T(G), g2 = T(2) * G;
  const size_t stride = size_t(W) + 1;
  for (i32 k = 0; k < n; ++k) {
    const i32 v0 = b + k, v1 = b + (k + 1 == n ? 0 : k + 1);
    const T ay = T(rs.y(v0)) * sy, by = T(rs.y(v1)) * sy;
    if (ay == by) continue;
    const T ax = T(rs.x(v0)) * sx, bx = T(rs.x(v1)) * sx;
    const i32 s = by > ay ? 1 : -1;
    T x0, y0, x1, y1;
    if (ay < by) { x0 = ax; y0 = ay; x1 = bx; y1 = by; }
    else         { x0 = bx; y0 = by; x1 = ax; y1 = ay; }
    const T dx = x1 - x0, dy = y1 - y0;  // dy > 0
    // rows r with y0 <= (2r+1) G < y1
    T rlo = -floordiv<T>(g - y0, g2), rhi = floordiv<T>(y1 - 1 - g, g2);
    if (rlo < 0) rlo = 0;
    if (rhi > T(H - 1)) rhi = T(H - 1);
    if (rlo > rhi) continue;
    const T gdy = g * dy;
    const i32 r0 = i32(rlo), r1 = i32(rhi);
    rmin = std::min(rmin, r0);
    rmax = std::max(rmax, r1);
    for (i32 r = r0; r <= r1; ++r) {
      const T py = T(2 * r + 1) * g;
      const T num = x0 * dy + (py - y0) * dx;  // crossing_x * dy
      // pixel c counts iff (2c+1) G dy < num ; q = max k with k G dy < num
      const T q = floordiv<T>(num - 1, gdy);
      const T cnt = floordiv<T>(q - 1, T(2)) + 1;  // number of counting pixels
      const i32 col = cnt < 0 ? 0 : (cnt > T(W) ? W : i32(cnt));
      diff[size_t(r) * stride + size_t(col)] -= s;
      cmin = std::min(cmin, col);
      cmax = std::max(cmax, col);
    }
  }
}

// Masks / statistics from the accumulated winding over the union region; restores the
// all-zero invariant of the scratch winding buffer. Specialised so the inner loop is
// branch-free.
template <bool MASKS, bool GT>
inline void final_pass(i32* wind, i32 W, i32 R0, i32 R1, i32 C0, i32 C1, bool own_wind,
                       const RasterOut& out, const u8* gt, RasterStats& st) {
  for (i32 r = R0; r <= R1; ++r) {
    i32* w = wind + size_t(r) * W;
    const size_t base = size_t(r) * W;
    std::uint32_t aeo = 0, anz = 0, multi = 0, ieo = 0, inz = 0;
    for (i32 c = C0; c < C1; ++c) {
      const i32 v = w[c];
      const u8 eo = u8(v & 1), nz = u8(v != 0);
      aeo += eo; anz += nz; multi += u8(unsigned(v + 1) > 2u);  // |v| >= 2
      if constexpr (MASKS) {
        if (out.evenodd) out.evenodd[base + c] = eo;
        if (out.nonzero) out.nonzero[base + c] = nz;
      }
      if constexpr (GT) {
        const u8 g = gt[base + c] != 0;
        ieo += eo & g; inz += nz & g;
      }
    }
    if (own_wind) std::fill(w + C0, w + C1, 0);
    st.area_evenodd += aeo; st.area_nonzero += anz; st.px_multi += multi;
    st.inter_evenodd += ieo; st.inter_nonzero += inz;
  }
}

template <class T>
inline void rasterize_impl(const Rings& rs, i32 W, i32 H, i64 G, RasterScratch& sc,
                           const RasterOut& out, const u8* gt, RasterStats& st) {
  sc.ensure(W, H);
  i32* diff = sc.diff.data();
  i32* wind = out.winding ? out.winding : sc.wind.data();
  const bool own_wind = out.winding == nullptr;
  if (!own_wind) std::fill(out.winding, out.winding + size_t(W) * H, 0);
  const size_t stride = size_t(W) + 1;
  i32 R0 = H, R1 = -1, C0 = W, C1 = 0;  // union region rows [R0,R1], cols [C0,C1)
  for (i32 ring = 0; ring < rs.nrings; ++ring) {
    i32 rmin = H, rmax = -1, cmin = W + 1, cmax = -1;
    scan_ring<T>(rs, ring, W, H, G, diff, rmin, rmax, cmin, cmax);
    i32 wmin = 0, wmax = 0;
    if (rmin <= rmax) {
      for (i32 r = rmin; r <= rmax; ++r) {
        i32* d = diff + size_t(r) * stride;
        i32* w = wind + size_t(r) * W;
        i32 acc = 0;
        for (i32 c = cmin; c < cmax; ++c) {
          acc += d[c];
          d[c] = 0;
          w[c] += acc;
          wmin = std::min(wmin, acc);
          wmax = std::max(wmax, acc);
        }
        d[cmax] = 0;
      }
      if (cmin < cmax) {
        R0 = std::min(R0, rmin); R1 = std::max(R1, rmax);
        C0 = std::min(C0, cmin); C1 = std::max(C1, cmax);
      }
    }
    if (out.ring_wmin) out.ring_wmin[ring] = wmin;
    if (out.ring_wmax) out.ring_wmax[ring] = wmax;
    st.rings_mixed_sign += (wmin < 0 && wmax > 0);
  }
  const size_t npx = size_t(W) * H;
  if (out.evenodd) std::memset(out.evenodd, 0, npx);
  if (out.nonzero) std::memset(out.nonzero, 0, npx);
  if (gt) {
    i64 a = 0;
    for (i32 r = 0; r < H; ++r) {
      const u8* g = gt + size_t(r) * W;
      std::uint32_t cnt = 0;
      for (i32 c = 0; c < W; ++c) cnt += g[c] != 0;
      a += cnt;
    }
    st.gt_area = a;
  }
  const bool masks = out.evenodd || out.nonzero;
  if (masks && gt) final_pass<true, true>(wind, W, R0, R1, C0, C1, own_wind, out, gt, st);
  else if (masks) final_pass<true, false>(wind, W, R0, R1, C0, C1, own_wind, out, gt, st);
  else if (gt) final_pass<false, true>(wind, W, R0, R1, C0, C1, own_wind, out, gt, st);
  else final_pass<false, false>(wind, W, R0, R1, C0, C1, own_wind, out, gt, st);
}

}  // namespace detail

// True iff int64 arithmetic is exact for this shape on this crop
// (all scaled coordinates and pixel centres below 2^30 in magnitude).
inline bool raster_fits_i64(const Rings& rs, i32 W, i32 H, i64 G) {
  const i64 dim = std::max(W, H);
  const i64 m = detail::max_abs_coord(rs);
  const i64 lim = i64(1) << 30;
  return m * 2 * dim < lim && 2 * dim * G < lim;
}

// Rasterize all rings together (fill rules apply to the union of rings).
// Preconditions: 1 <= W, H <= kMaxDim, 1 <= G <= kMaxGrid. gt (optional) is an
// H*W 0/1 mask used only to count intersections for IoU.
inline RasterStats rasterize(const Rings& rs, i32 W, i32 H, i64 G, RasterScratch& sc,
                             const RasterOut& out, const u8* gt = nullptr) {
  RasterStats st;
  if (raster_fits_i64(rs, W, H, G)) detail::rasterize_impl<i64>(rs, W, H, G, sc, out, gt, st);
  else detail::rasterize_impl<i128>(rs, W, H, G, sc, out, gt, st);
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
  kHoleOutsideOuter,      // hole with a witness point strictly outside its outer ring
  kNestedHoles,           // # sibling holes this hole overlaps / is nested with
  kNestedOuters,          // # other outer parts this outer overlaps / is nested with
  kWrongWinding,          // outer with area < 0 or hole with area > 0
  kOutOfRange,            // # vertices with x or y outside [0, G]
  kNumDefects
};

inline constexpr const char* kDefectNames[kNumDefects] = {
    "fewer_than_3_points", "zero_area",          "self_crossing",
    "vertex_pinch",        "collinear_overlap",  "spike",
    "ring_ring_crossing",  "ring_ring_touching", "hole_outside_outer",
    "nested_holes",        "nested_outers",      "wrong_winding",
    "out_of_range"};

enum class PairMethod { Sweep, Brute };

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
  i32 at(i32 r, Defect d) const { return counts[size_t(r) * kNumDefects + d]; }
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
                              ValidateScratch& s, PairMethod method) {
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
  if (method == PairMethod::Brute) {
    for (i32 i = 0; i < E; ++i)
      for (i32 j = i + 1; j < E; ++j)
        if (boxes_meet(ed[i], ed[j])) record_pair<T>(ed[i], ed[j], out);
  } else {  // sort-and-sweep on x, test y-overlap
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
// index) for a hole of ring p (the binding validates this). G is the grid size
// used for out_of_range. Never repairs anything.
inline void validate(const Rings& rs, const i32* roles, i64 G, ValidateResult& out,
                     ValidateScratch& s, PairMethod method = PairMethod::Sweep) {
  const i32 R = rs.nrings;
  out.reset(R);
  s.px.clear(); s.py.clear();
  s.start.assign(R, 0); s.len.assign(R, 0); s.fewer.assign(R, 0);
  s.bx0.assign(R, 0); s.bx1.assign(R, 0); s.by0.assign(R, 0); s.by1.assign(R, 0);
  i64 M = 0;
  for (i32 r = 0; r < R; ++r) {
    const i32 b = rs.begin(r), n = rs.size(r);
    const i32 st = i32(s.px.size());
    s.start[r] = st;
    i32 oor = 0;
    for (i32 k = 0; k < n; ++k) {
      const i64 x = rs.x(b + k), y = rs.y(b + k);
      oor += (x < 0 || x > G || y < 0 || y > G);
      M = std::max(M, std::max(detail::iabs(x), detail::iabs(y)));
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
  if (M < (i64(1) << 28)) detail::validate_geometry<i64>(rs, roles, out, s, method);
  else detail::validate_geometry<i128>(rs, roles, out, s, method);
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
inline void score_batch(const Rollout* in, i64 n, i64 G, int threads, RolloutScore* out) {
  if (threads <= 1 || n <= 1) {
    Workspace ws;
    for (i64 i = 0; i < n; ++i) score_one(in[i], G, ws, out[i]);
    return;
  }
  std::atomic<i64> next{0};
  auto work = [&]() {
    Workspace ws;
    for (i64 i; (i = next.fetch_add(1, std::memory_order_relaxed)) < n;) score_one(in[i], G, ws, out[i]);
  };
  const int nt = int(std::min<i64>(threads, n));
  std::vector<std::thread> pool;
  pool.reserve(nt - 1);
  for (int t = 1; t < nt; ++t) pool.emplace_back(work);
  work();
  for (auto& th : pool) th.join();
}

}  // namespace polycore
