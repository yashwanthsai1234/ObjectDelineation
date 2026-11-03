// gpucheck.cpp -- the single-source checker of gpucheck.h. Every algorithm is written once, as __host__ __device__
// functors run by for_each over independent items (rows, rings, edges, shapes); for_each is a plain loop when the
// C++ compiler builds this file and a kernel launch when nvcc does. The only other CPU/GPU differences are memory
// (Buf) and four atomic helpers, all in the "platform" section below.
//
// Polygons (port of cpp/polycore/polycore.hpp, same exact integer arithmetic, int64 or __int128 chosen per shape
// by the same bounds):
//  * raster: one item per pixel row. Each edge crossing the row's centre line adds its sign to the winding of the
//    centres left of the crossing; a suffix sum over a 256-column window gives the winding, so no sort and no
//    scratch memory. Even-odd and nonzero come from the same winding (one pass).
//  * validate: one item per ring (clean the ring, count repeats/area/range), one per edge (classify it against
//    every later edge of the shape; polycore's sweep only prunes pairs whose boxes do not meet, which classify
//    rejects anyway), one per shape (ring-ring counts, winding, containment by exact witness points).
// Masks (port of objdel.checker.polygonize + bitquad, the array-only design):
//  * every unit boundary edge is a slot (vertex, direction); its successor is read from the 2x2 bit-quad at its end
//    (left turn at saddles: 8-connected foreground, 4-connected background);
//  * ring label = minimum slot on the cycle (pointer jumping); the minimum slot is polygonize's canonical start;
//  * corner index = number of corners from the start (weighted list ranking on the predecessor);
//  * rings are ordered as polygonize returns them (by their first east edge in raster order) and their corners
//    scattered to their output positions.
//
// References: the headers of cpp/polycore/polycore.hpp and src/objdel/checker/polygonize.py; CUDA __int128 in device
// code (CUDA >= 11.5): https://docs.nvidia.com/cuda/cuda-c-programming-guide/#supported-types ; unified memory:
// https://docs.nvidia.com/cuda/cuda-c-programming-guide/#unified-memory-programming
#include "gpucheck.h"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string>

#include "polycore/polycore.hpp"  // Defect, kNumDefects, kMaxRings (constants only; the kernels are below)

#ifdef __CUDACC__
#include <cuda_runtime.h>
#define HD __host__ __device__
#else
#define HD
#endif

namespace gpucheck {
namespace {

using i128 = __int128;
using u64 = std::uint64_t;
using polycore::kMaxRings;
using polycore::kNumDefects;
constexpr i32 kPairs = kMaxRings * kMaxRings;

// ===================================================================== platform
#ifdef __CUDACC__
void cuda_check(cudaError_t e) {
  if (e != cudaSuccess) throw std::runtime_error(std::string("CUDA: ") + cudaGetErrorString(e));
}

template <class F>
__global__ void kernel(i64 n, F f) {
  const i64 i = i64(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < n) f(i);
}

template <class F>
void for_each(i64 n, const F& f) {
  constexpr int kThreads = 128;
  if (n <= 0) return;
  kernel<<<unsigned((n + kThreads - 1) / kThreads), kThreads>>>(n, f);
  cuda_check(cudaGetLastError());
}

void sync() { cuda_check(cudaDeviceSynchronize()); }
#else
template <class F>
void for_each(i64 n, const F& f) {
  for (i64 i = 0; i < n; ++i) f(i);
}

void sync() {}
#endif

// Zero-filled array readable by the host and the kernels (unified memory on the GPU). Host access only after sync().
template <class T>
struct Buf {
  T* p = nullptr;
  explicit Buf(i64 n, const T* src = nullptr) {
    const size_t bytes = sizeof(T) * size_t(n > 0 ? n : 1);
#ifdef __CUDACC__
    cuda_check(cudaMallocManaged(&p, bytes));
    cuda_check(cudaMemset(p, 0, bytes));
    sync();
#else
    p = static_cast<T*>(::operator new(bytes));
    std::memset(p, 0, bytes);
#endif
    if (src && n > 0) std::memcpy(p, src, sizeof(T) * size_t(n));
  }
  ~Buf() {
#ifdef __CUDACC__
    cudaFree(p);
#else
    ::operator delete(p);
#endif
  }
  Buf(const Buf&) = delete;
  Buf& operator=(const Buf&) = delete;
};

HD inline void atomic_add(i32* p, i32 v) {
#ifdef __CUDA_ARCH__
  atomicAdd(p, v);
#else
  *p += v;
#endif
}

HD inline void atomic_add(i64* p, i64 v) {  // two's complement: unsigned addition is signed addition
#ifdef __CUDA_ARCH__
  atomicAdd(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(v));
#else
  *p += v;
#endif
}

HD inline i32 atomic_inc(i32* p) {  // returns the old value
#ifdef __CUDA_ARCH__
  return atomicAdd(p, 1);
#else
  return (*p)++;
#endif
}

HD inline void atomic_min(i32* p, i32 v) {
#ifdef __CUDA_ARCH__
  atomicMin(p, v);
#else
  if (v < *p) *p = v;
#endif
}

// ===================================================================== helpers
template <class T>
HD inline T mn(T a, T b) { return b < a ? b : a; }
template <class T>
HD inline T mx(T a, T b) { return a < b ? b : a; }

template <class T>
HD inline T floordiv(T a, T b) {  // b != 0
  T q = a / b;
  if ((a % b != 0) && ((a < 0) != (b < 0))) --q;
  return q;
}

// Segment of item i in a prefix-sum array start[0..n]: the largest s < n with start[s] <= i.
HD inline i64 segment(const i64* start, i64 n, i64 i) {
  i64 lo = 0, hi = n;
  while (hi - lo > 1) {
    const i64 mid = (lo + hi) / 2;
    if (start[mid] <= i) lo = mid;
    else hi = mid;
  }
  return lo;
}

HD inline i64 iabs(i64 v) { return v < 0 ? -v : v; }

// ===================================================================== polygons
struct Shapes {  // device view of Polygons plus derived offsets
  const i32* xy;
  const i64 *ring_start, *shape_start, *row_start, *mask_start;
  const i32 *roles, *W, *H;
  i64 n, G;
};

// polycore's precision choices, per shape: rasterize (raster_fits_i64) and validate (max |coord| < 2^28).
struct Precision {
  Shapes s;
  u8 *raster_wide, *validate_wide;
  HD void operator()(i64 sh) const {
    i64 m = 0;
    for (i64 v = s.ring_start[s.shape_start[sh]]; v < s.ring_start[s.shape_start[sh + 1]]; ++v)
      m = mx(m, mx(iabs(s.xy[2 * v]), iabs(s.xy[2 * v + 1])));
    const i64 dim = mx<i64>(s.W[sh], s.H[sh]), lim = i64(1) << 30;
    raster_wide[sh] = !(m * 2 * dim < lim && 2 * dim * s.G < lim);
    validate_wide[sh] = !(m < (i64(1) << 28));
  }
};

// ---------------------------------------------------------------- raster
constexpr i32 kChunk = 256;  // columns per suffix-sum window (local memory per thread: 1 KiB)

template <class T>
HD void raster_row(const Shapes& s, i64 sh, i32 row, u8* eo, u8* nz) {
  const i32 W = s.W[sh];
  const T sx = T(2) * W, sy = T(2) * s.H[sh], g = T(s.G), yc = T(2 * row + 1) * g;
  for (i32 c0 = 0; c0 < W; c0 += kChunk) {
    const i32 c1 = mn(W, c0 + kChunk);
    i32 buf[kChunk];
    for (i32 k = 0; k < c1 - c0; ++k) buf[k] = 0;
    i32 right = 0;  // signs of crossings right of the whole window
    for (i64 ring = s.shape_start[sh]; ring < s.shape_start[sh + 1]; ++ring) {
      const i64 b = s.ring_start[ring], n = s.ring_start[ring + 1] - b;
      for (i64 k = 0; k < n; ++k) {
        const i64 v0 = b + k, v1 = b + (k + 1 == n ? 0 : k + 1);
        const T ay = T(s.xy[2 * v0 + 1]) * sy, by = T(s.xy[2 * v1 + 1]) * sy;
        if (ay == by) continue;  // horizontal edges never count
        const bool down = by > ay;  // edge goes down the screen: +1
        const T y0 = down ? ay : by, dy = down ? by - ay : ay - by;
        if (yc < y0 || yc >= y0 + dy) continue;  // half-open: y0 <= yc < y0 + dy
        const T ax = T(s.xy[2 * v0]) * sx, bx = T(s.xy[2 * v1]) * sx;
        const T x0 = down ? ax : bx, dx = (down ? bx : ax) - x0;
        // col = number of centres strictly left of the crossing, clamped to [0, W] (polycore collect_keys)
        const T q = floordiv<T>(x0 * dy + (yc - y0) * dx - 1, g * dy);
        const T cnt = floordiv<T>(q - 1, T(2)) + 1;
        const i32 col = cnt < 0 ? 0 : (cnt > T(W) ? W : i32(cnt));
        const i32 sign = down ? 1 : -1;
        if (col >= c1) right += sign;
        else if (col > c0) buf[col - 1 - c0] += sign;  // adds to centres c0 .. col - 1
      }
    }
    i32 w = right;  // winding(c) = sum of the signs of crossings right of centre c
    for (i32 c = c1 - 1; c >= c0; --c) {
      w += buf[c - c0];
      eo[c] = u8(w & 1);
      nz[c] = u8(w != 0);
    }
  }
}

struct RasterRows {
  Shapes s;
  const u8* raster_wide;
  u8 *eo, *nz;
  HD void operator()(i64 i) const {
    const i64 sh = segment(s.row_start, s.n, i);
    const i32 row = i32(i - s.row_start[sh]);
    const i64 at = s.mask_start[sh] + i64(row) * s.W[sh];
    if (raster_wide[sh]) raster_row<i128>(s, sh, row, eo + at, nz + at);
    else raster_row<i64>(s, sh, row, eo + at, nz + at);
  }
};

// ---------------------------------------------------------------- validate
struct Rings {  // cleaned rings (consecutive and closing repeats dropped), stored at their input offsets
  i32* cxy;     // clean vertices
  i32* len;     // clean vertex count per ring
  u8* fewer;    // < 3 distinct vertices
  i64* area2;   // saturated twice signed area of the clean ring
  i32* box;     // x0, x1, y0, y1 per ring
  i32* counts;  // (rings, 13) per-ring defect counts
  i32* pairs;   // per shape: cross (10 x 10), touch (10 x 10) edge-pair counts between rings
};

HD inline i64 saturate(i128 v) {
  const i128 hi = i128(0x7FFFFFFFFFFFFFFFll), lo = -hi - 1;
  return v > hi ? i64(hi) : (v < lo ? i64(lo) : i64(v));
}

struct CleanRings {
  Shapes s;
  Rings r;
  HD void operator()(i64 ring) const {
    const i64 b = s.ring_start[ring], n = s.ring_start[ring + 1] - b;
    const i32* p = s.xy + 2 * b;
    i32* c = r.cxy + 2 * b;
    i32* cnt = r.counts + ring * kNumDefects;
    i32 m = 0, oor = 0;
    for (i64 k = 0; k < n; ++k) {
      const i64 x = p[2 * k], y = p[2 * k + 1];
      oor += (x < 0 || x > s.G || y < 0 || y > s.G);
      if (k == 0 || x != p[2 * k - 2] || y != p[2 * k - 1]) c[2 * m] = i32(x), c[2 * m + 1] = i32(y), ++m;
    }
    if (m > 1 && c[0] == c[2 * m - 2] && c[1] == c[2 * m - 1]) --m;
    r.len[ring] = m;
    cnt[polycore::kOutOfRange] = oor;
    i32* box = r.box + 4 * ring;
    for (i32 k = 0; k < m; ++k) {
      box[0] = k ? mn(box[0], c[2 * k]) : c[0];
      box[1] = k ? mx(box[1], c[2 * k]) : c[0];
      box[2] = k ? mn(box[2], c[2 * k + 1]) : c[1];
      box[3] = k ? mx(box[3], c[2 * k + 1]) : c[1];
    }
    i32 distinct = m > 0 ? 1 : 0, second = -1;
    for (i32 k = 1; k < m && distinct < 3; ++k) {
      if (c[2 * k] == c[0] && c[2 * k + 1] == c[1]) continue;
      if (second < 0) { second = k; distinct = 2; continue; }
      if (c[2 * k] == c[2 * second] && c[2 * k + 1] == c[2 * second + 1]) continue;
      distinct = 3;
    }
    i128 a2 = 0;
    for (i32 k = 0; k < m; ++k) {
      const i32 k1 = k + 1 == m ? 0 : k + 1;
      a2 += i128(c[2 * k]) * c[2 * k1 + 1] - i128(c[2 * k1]) * c[2 * k + 1];
    }
    r.area2[ring] = saturate(a2);
    r.fewer[ring] = distinct < 3;
    if (distinct < 3) cnt[polycore::kFewerThan3Points] = 1;
    else if (a2 == 0) cnt[polycore::kZeroArea] = 1;
  }
};

struct Edge {
  i64 ax, ay, bx, by, xmin, xmax, ymin, ymax;
  i32 ring, idx, n;  // ring local to the shape, index on the clean ring, clean ring length
};

HD inline Edge make_edge(const Rings& r, const Shapes& s, i64 ring, i32 local, i32 k) {
  const i32* c = r.cxy + 2 * s.ring_start[ring];
  const i32 n = r.len[ring], k1 = k + 1 == n ? 0 : k + 1;
  Edge e;
  e.ax = c[2 * k], e.ay = c[2 * k + 1], e.bx = c[2 * k1], e.by = c[2 * k1 + 1];
  e.xmin = mn(e.ax, e.bx), e.xmax = mx(e.ax, e.bx), e.ymin = mn(e.ay, e.by), e.ymax = mx(e.ay, e.by);
  e.ring = local, e.idx = k, e.n = n;
  return e;
}

template <class T>
HD inline int orient(i64 ax, i64 ay, i64 bx, i64 by, i64 cx, i64 cy) {
  const T v = T(bx - ax) * T(cy - ay) - T(by - ay) * T(cx - ax);
  return (v > 0) - (v < 0);
}

HD inline bool in_box(i64 ax, i64 ay, i64 bx, i64 by, i64 px, i64 py) {
  return mn(ax, bx) <= px && px <= mx(ax, bx) && mn(ay, by) <= py && py <= mx(ay, by);
}

enum PairKind { kNone = 0, kCross, kPinch, kOverlap, kSpikeKind };

// polycore detail::classify, unchanged.
template <class T>
HD inline PairKind classify(const Edge& e, const Edge& f) {
  if (e.ring == f.ring) {
    const i32 n = e.n;
    const bool e_then_f = (f.idx == (e.idx + 1 == n ? 0 : e.idx + 1));
    const bool f_then_e = (e.idx == (f.idx + 1 == n ? 0 : f.idx + 1));
    if (e_then_f || f_then_e) {
      const Edge& a = e_then_f ? e : f;  // edge ending at the shared vertex S
      const Edge& b = e_then_f ? f : e;  // edge starting at S
      const i64 sx = a.bx, sy = a.by, px = a.ax, py = a.ay, qx = b.bx, qy = b.by;
      if (orient<T>(sx, sy, px, py, qx, qy) == 0 && (in_box(sx, sy, px, py, qx, qy) || in_box(sx, sy, qx, qy, px, py)))
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
    const i64 lo = use_x ? mx(e.xmin, f.xmin) : mx(e.ymin, f.ymin);
    const i64 hi = use_x ? mn(e.xmax, f.xmax) : mn(e.ymax, f.ymax);
    if (hi > lo) return kOverlap;
  }
  return kPinch;
}

// One item per input vertex slot: the clean edge starting there (if any) against every later edge of its shape.
struct EdgePairs {
  Shapes s;
  Rings r;
  const u8* validate_wide;
  template <class T>
  HD void run(i64 sh, i64 ring, i32 k) const {
    const i64 r0 = s.shape_start[sh], r1 = s.shape_start[sh + 1];
    const Edge e = make_edge(r, s, ring, i32(ring - r0), k);
    i32* pc = r.pairs + sh * 2 * kPairs;
    for (i64 q = ring; q < r1; ++q) {
      if (r.fewer[q]) continue;
      for (i32 j = q == ring ? k + 1 : 0; j < r.len[q]; ++j) {
        const Edge f = make_edge(r, s, q, i32(q - r0), j);
        if (!(e.xmin <= f.xmax && f.xmin <= e.xmax && e.ymin <= f.ymax && f.ymin <= e.ymax)) continue;
        const PairKind kind = classify<T>(e, f);
        if (kind == kNone) continue;
        if (q == ring) {
          const int d = kind == kCross ? polycore::kSelfCrossing
                      : kind == kPinch ? polycore::kVertexPinch
                      : kind == kOverlap ? polycore::kCollinearOverlap : polycore::kSpike;
          atomic_add(r.counts + ring * kNumDefects + d, 1);
        } else {
          i32* m = pc + (kind == kCross ? 0 : kPairs);
          atomic_add(m + e.ring * kMaxRings + f.ring, 1);
          atomic_add(m + f.ring * kMaxRings + e.ring, 1);
        }
      }
    }
  }
  HD void operator()(i64 v) const {
    const i64 nr = s.shape_start[s.n];
    const i64 ring = segment(s.ring_start, nr, v);
    const i32 k = i32(v - s.ring_start[ring]);
    if (r.fewer[ring] || k >= r.len[ring]) return;
    const i64 sh = segment(s.shape_start, s.n, ring);
    if (validate_wide[sh]) run<i128>(sh, ring, k);
    else run<i64>(sh, ring, k);
  }
};

// Ring-ring counts, winding and containment of one shape (polycore detail::validate_geometry, second half).
struct ShapeTotals {
  Shapes s;
  Rings r;
  const u8* validate_wide;
  i64* defects;

  // Point in doubled coordinates against clean ring g: -1 strictly outside, 0 on the boundary, 1 inside.
  template <class T>
  HD int locate2(i64 g, i64 qx, i64 qy) const {
    const i32* box = r.box + 4 * g;
    if (qx < 2 * i64(box[0]) || qx > 2 * i64(box[1]) || qy < 2 * i64(box[2]) || qy > 2 * i64(box[3])) return -1;
    const i32* c = r.cxy + 2 * s.ring_start[g];
    const i32 n = r.len[g];
    bool in = false;
    for (i32 k = 0; k < n; ++k) {
      const i32 k1 = k + 1 == n ? 0 : k + 1;
      const i64 ax = 2 * i64(c[2 * k]), ay = 2 * i64(c[2 * k + 1]), bx = 2 * i64(c[2 * k1]), by = 2 * i64(c[2 * k1 + 1]);
      if (in_box(ax, ay, bx, by, qx, qy) && orient<T>(ax, ay, bx, by, qx, qy) == 0) return 0;
      if ((ay <= qy && qy < by) || (by <= qy && qy < ay)) {
        const T num = T(ax - qx) * T(by - ay) + T(qy - ay) * T(bx - ax);
        if (by > ay ? num > 0 : num < 0) in = !in;
      }
    }
    return in ? 1 : -1;
  }

  // Region of outer ring o (local) with its holes, even-odd over all of them.
  template <class T>
  HD int locate2_region(i64 r0, i32 R, i32 o, i64 qx, i64 qy) const {
    const int loc = locate2<T>(r0 + o, qx, qy);
    if (loc == 0) return 0;
    bool in = loc > 0;
    for (i32 h = 0; h < R; ++h) {
      if (s.roles[r0 + h] != o || r.fewer[r0 + h]) continue;
      const int lh = locate2<T>(r0 + h, qx, qy);
      if (lh == 0) return 0;
      if (lh > 0) in = !in;
    }
    return in ? 1 : -1;
  }

  // Does a witness of ring a (vertex or edge midpoint) locate to `want` against ring b (or b's region)?
  // single: a neither crosses nor touches that boundary, so its first vertex decides.
  template <class T>
  HD bool any_witness(i64 r0, i32 R, i32 a, int want, bool single, i32 b, bool region) const {
    const i32* c = r.cxy + 2 * s.ring_start[r0 + a];
    const i32 n = r.len[r0 + a];
    auto loc = [&](i64 x, i64 y) { return region ? locate2_region<T>(r0, R, b, x, y) : locate2<T>(r0 + b, x, y); };
    if (single) return n > 0 && loc(2 * i64(c[0]), 2 * i64(c[1])) == want;
    for (i32 k = 0; k < n; ++k) {
      const i32 k1 = k + 1 == n ? 0 : k + 1;
      if (loc(2 * i64(c[2 * k]), 2 * i64(c[2 * k + 1])) == want) return true;
      if (loc(i64(c[2 * k]) + c[2 * k1], i64(c[2 * k + 1]) + c[2 * k1 + 1]) == want) return true;
    }
    return false;
  }

  template <class T>
  HD void run(i64 sh) const {
    const i64 r0 = s.shape_start[sh];
    const i32 R = i32(s.shape_start[sh + 1] - r0);
    const i32* roles = s.roles + r0;
    const i32 *pcross = r.pairs + sh * 2 * kPairs, *ptouch = pcross + kPairs;
    auto at = [&](i32 ring, int d) -> i32& { return r.counts[(r0 + ring) * kNumDefects + d]; };
    auto fewer = [&](i32 ring) { return r.fewer[r0 + ring] != 0; };
    auto apart = [&](i32 a, i32 b) { return pcross[a * kMaxRings + b] == 0 && ptouch[a * kMaxRings + b] == 0; };
    auto apart_region = [&](i32 a, i32 o) {  // a vs outer o and all of o's holes
      if (!apart(a, o)) return false;
      for (i32 h = 0; h < R; ++h)
        if (roles[h] == o && !apart(a, h)) return false;
      return true;
    };
    i64* tot = defects + sh * kNumDefects;
    for (i32 a = 0; a < R; ++a)
      for (i32 b = 0; b < R; ++b) {
        at(a, polycore::kRingRingCrossing) += pcross[a * kMaxRings + b];
        at(a, polycore::kRingRingTouching) += ptouch[a * kMaxRings + b];
      }
    for (i32 a = 0; a < R; ++a) {
      if (fewer(a)) continue;
      const i64 a2 = r.area2[r0 + a];
      if ((roles[a] < 0 && a2 < 0) || (roles[a] >= 0 && a2 > 0)) at(a, polycore::kWrongWinding) = 1;
    }
    for (i32 h = 0; h < R; ++h) {
      const i32 p = roles[h];
      if (p < 0 || fewer(h) || fewer(p)) continue;
      if (any_witness<T>(r0, R, h, -1, apart(h, p), p, false)) at(h, polycore::kHoleOutsideOuter) = 1;
    }
    for (i32 a = 0; a < R; ++a) {
      if (fewer(a)) continue;
      for (i32 b = a + 1; b < R; ++b) {
        if (fewer(b)) continue;
        if (roles[a] >= 0 && roles[a] == roles[b]) {  // sibling holes
          const bool single = apart(a, b);
          if (any_witness<T>(r0, R, b, 1, single, a, false) || any_witness<T>(r0, R, a, 1, single, b, false)) {
            at(a, polycore::kNestedHoles) += 1, at(b, polycore::kNestedHoles) += 1;
            tot[polycore::kNestedHoles] += 1;
          }
        } else if (roles[a] < 0 && roles[b] < 0) {  // two outers
          if (any_witness<T>(r0, R, b, 1, apart_region(b, a), a, true) ||
              any_witness<T>(r0, R, a, 1, apart_region(a, b), b, true)) {
            at(a, polycore::kNestedOuters) += 1, at(b, polycore::kNestedOuters) += 1;
            tot[polycore::kNestedOuters] += 1;
          }
        }
      }
    }
    for (i32 a = 0; a < R; ++a)
      for (int d = 0; d < kNumDefects; ++d)
        if (d != polycore::kRingRingCrossing && d != polycore::kRingRingTouching && d != polycore::kNestedHoles &&
            d != polycore::kNestedOuters)
          tot[d] += at(a, d);
    for (i32 a = 0; a < R; ++a)
      for (i32 b = a + 1; b < R; ++b) {
        tot[polycore::kRingRingCrossing] += pcross[a * kMaxRings + b];
        tot[polycore::kRingRingTouching] += ptouch[a * kMaxRings + b];
      }
  }

  HD void operator()(i64 sh) const {
    if (validate_wide[sh]) run<i128>(sh);
    else run<i64>(sh);
  }
};

// ===================================================================== masks
// Slot of the unit edge leaving vertex (x, y) in direction d (0 east, 1 south, 2 west, 3 north), local to a mask:
// (y (W + 1) + x) 4 + d, which is also polygonize's edge key, so the minimum slot of a ring is its canonical start.
// The bit-quad at vertex (x, y) has TL = pixel (x-1, y-1) = 1, TR = (x, y-1) = 2, BL = (x-1, y) = 4, BR = (x, y) = 8,
// as in checker.bitquad. Foreground lies right of every edge (clockwise around parts, counter-clockwise in holes).
struct MaskView {
  const u8* bits;
  const i64 *start, *slot_start, *vrow_start;
  const i32 *W, *H;
  i64 n;

  HD bool px(i64 m, i32 x, i32 y) const {
    return x >= 0 && y >= 0 && x < W[m] && y < H[m] && bits[start[m] + i64(y) * W[m] + x] != 0;
  }
  HD int quad(i64 m, i32 x, i32 y) const {
    return px(m, x - 1, y - 1) | px(m, x, y - 1) << 1 | px(m, x - 1, y) << 2 | px(m, x, y) << 3;
  }
};

HD inline bool has_edge(int code, int d) {
  switch (d) {
    case 0: return (code & 8) && !(code & 2);  // top side of BR
    case 1: return (code & 4) && !(code & 8);  // right side of BL
    case 2: return (code & 1) && !(code & 4);  // bottom side of TL
    default: return (code & 2) && !(code & 1);  // left side of TR
  }
}

HD inline i32 step_x(int d) { return d == 0 ? 1 : (d == 2 ? -1 : 0); }
HD inline i32 step_y(int d) { return d == 1 ? 1 : (d == 3 ? -1 : 0); }
HD inline bool is_corner(const i32* pred, i32 e) { return (e & 3) != (pred[e] & 3); }

struct Slots {     // dense per-slot arrays of one chunk of masks, indexed slot_start[m] + local slot
  i32* succ;       // successor slot, -1 = no edge
  i32* pred;
  i32 *lab[2], *jmp[2], *w[2];  // ping-pong buffers of the pointer-jumping rounds
  i32* rcount;     // per ring (at its start slot): corners
  i32* rkey;       // per ring: first east edge in raster order (polygonize's ring order)
  i64* rarea;      // per ring: twice the signed area
  i64* rbase;      // per ring: output offset of its first corner
};

// The outer loops below walk one vertex row (mask m, y) per item; f(m, base, x, y, local slot e).
template <class F>
HD inline void for_row_slots(const MaskView& v, i64 i, const Slots& sl, F f) {
  const i64 m = segment(v.vrow_start, v.n, i);
  const i32 y = i32(i - v.vrow_start[m]), W = v.W[m];
  const i64 base = v.slot_start[m];
  for (i32 x = 0; x <= W; ++x)
    for (int d = 0; d < 4; ++d) {
      const i32 e = (y * (W + 1) + x) * 4 + d;
      if (sl.succ[base + e] >= 0) f(m, base, x, y, e);
    }
}

struct InitSlots {
  MaskView v;
  Slots sl;
  i64* quads;
  HD void operator()(i64 i) const {
    const i64 m = segment(v.vrow_start, v.n, i);
    const i32 y = i32(i - v.vrow_start[m]), W = v.W[m];
    const i64 base = v.slot_start[m];
    i64 n1 = 0, n2 = 0, nd = 0, n3 = 0, n4 = 0;
    for (i32 x = 0; x <= W; ++x) {
      const int code = v.quad(m, x, y);
      const int bits = (code & 1) + (code >> 1 & 1) + (code >> 2 & 1) + (code >> 3);
      if (code == 6 || code == 9) ++nd;
      else if (bits == 1) ++n1;
      else if (bits == 2) ++n2;
      else if (bits == 3) ++n3;
      else if (bits == 4) ++n4;
      for (int d = 0; d < 4; ++d) {
        const i32 e = (y * (W + 1) + x) * 4 + d;
        i32 next = -1;
        if (has_edge(code, d)) {  // successor: the left turn if it exists (saddles), else the only edge
          const i32 x1 = x + step_x(d), y1 = y + step_y(d);
          const int c1 = v.quad(m, x1, y1);
          const int turns[3] = {(d + 3) & 3, d, (d + 1) & 3};
          for (int t : turns)
            if (next < 0 && has_edge(c1, t)) next = (y1 * (W + 1) + x1) * 4 + t;
        }
        sl.succ[base + e] = next;
        sl.lab[0][base + e] = e;
        sl.jmp[0][base + e] = next;
        sl.rkey[base + e] = 0x7FFFFFFF;
      }
    }
    i64* q = quads + 5 * m;
    atomic_add(q + 0, n1), atomic_add(q + 1, n2), atomic_add(q + 2, nd), atomic_add(q + 3, n3), atomic_add(q + 4, n4);
  }
};

// Label round: lab = min(lab, lab[jmp]); jmp = jmp[jmp]. After k rounds lab covers 2^k slots of the cycle.
struct LabelRound {
  MaskView v;
  Slots sl;
  int in;
  HD void operator()(i64 i) const {
    const i32 *lab = sl.lab[in], *jmp = sl.jmp[in];
    i32 *lab2 = sl.lab[1 - in], *jmp2 = sl.jmp[1 - in];
    for_row_slots(v, i, sl, [&](i64, i64 base, i32, i32, i32 e) {
      const i32 j = jmp[base + e];
      lab2[base + e] = mn(lab[base + e], lab[base + j]);
      jmp2[base + e] = jmp[base + j];
    });
  }
};

struct Pred {
  MaskView v;
  Slots sl;
  HD void operator()(i64 i) const {
    for_row_slots(v, i, sl, [&](i64, i64 base, i32, i32, i32 e) { sl.pred[base + sl.succ[base + e]] = e; });
  }
};

// Ranking on pred, stopped at the ring's start: w(e) = corners from the start up to e's predecessor = e's corner index.
struct RankInit {
  MaskView v;
  Slots sl;
  const i32* lab;
  HD void operator()(i64 i) const {
    for_row_slots(v, i, sl, [&](i64, i64 base, i32, i32, i32 e) {
      const bool start = lab[base + e] == e;
      const i32 p = sl.pred[base + e];
      sl.jmp[0][base + e] = start ? e : p;
      sl.w[0][base + e] = start ? 0 : is_corner(sl.pred + base, p);
    });
  }
};

struct RankRound {
  MaskView v;
  Slots sl;
  int in;
  HD void operator()(i64 i) const {
    const i32 *w = sl.w[in], *jmp = sl.jmp[in];
    i32 *w2 = sl.w[1 - in], *jmp2 = sl.jmp[1 - in];
    for_row_slots(v, i, sl, [&](i64, i64 base, i32, i32, i32 e) {
      const i32 j = jmp[base + e];
      w2[base + e] = w[base + e] + w[base + j];
      jmp2[base + e] = jmp[base + j];
    });
  }
};

struct RingStats {
  MaskView v;
  Slots sl;
  const i32* lab;
  i32 *nrings, *rec;  // rec: (mask chunk index, start slot) of every ring, in any order
  HD void operator()(i64 i) const {
    for_row_slots(v, i, sl, [&](i64 m, i64 base, i32 x, i32 y, i32 e) {
      const i32 root = lab[base + e], d = e & 3, x1 = x + step_x(d), y1 = y + step_y(d);
      atomic_add(sl.rarea + base + root, i64(x) * y1 - i64(x1) * y);
      if (is_corner(sl.pred + base, e)) atomic_add(sl.rcount + base + root, 1);
      if (d == 0) atomic_min(sl.rkey + base + root, y * (v.W[m] + 1) + x);
      if (root == e) {
        const i32 k = atomic_inc(nrings);
        rec[2 * k] = i32(m), rec[2 * k + 1] = e;
      }
    });
  }
};

struct ScatterCorners {
  MaskView v;
  Slots sl;
  const i32 *lab, *w;
  i32* pts;
  HD void operator()(i64 i) const {
    for_row_slots(v, i, sl, [&](i64, i64 base, i32 x, i32 y, i32 e) {
      if (!is_corner(sl.pred + base, e)) return;
      const i64 k = sl.rbase[base + lab[base + e]] + w[base + e];
      pts[2 * k] = x, pts[2 * k + 1] = y;
    });
  }
};

int jump_rounds(i64 n) {  // polygonize._jumps: ceil(log2(n)) + 1 rounds cover any cycle of <= n edges
  int k = 0;
  while ((i64(1) << k) < n) ++k;
  return k + 1;
}

constexpr i64 kChunkSlots = i64(1) << 23;  // slots per batch of masks (~0.5 GB of scratch)

void check_mask_chunk(const Masks& in, i64 m0, i64 m1, MaskResult& out) {
  const i64 n = m1 - m0;
  std::vector<i64> start(n + 1), slot_start(n + 1, 0), vrow_start(n + 1, 0);
  i64 max_slots = 1;
  for (i64 k = 0; k < n; ++k) {
    const i64 W = in.W[m0 + k], H = in.H[m0 + k];
    start[k] = in.start[m0 + k] - in.start[m0];
    slot_start[k + 1] = slot_start[k] + 4 * (W + 1) * (H + 1);
    vrow_start[k + 1] = vrow_start[k] + H + 1;
    max_slots = std::max(max_slots, 4 * (W + 1) * (H + 1));
  }
  start[n] = in.start[m1] - in.start[m0];
  const i64 S = slot_start[n], rows = vrow_start[n];
  Buf<u8> bits(start[n], in.bits + in.start[m0]);
  Buf<i64> d_start(n + 1, start.data()), d_slot(n + 1, slot_start.data()), d_vrow(n + 1, vrow_start.data());
  Buf<i32> W(n, in.W + m0), H(n, in.H + m0);
  Buf<i32> succ(S), pred(S), lab0(S), lab1(S), jmp0(S), jmp1(S), w0(S), w1(S), rcount(S), rkey(S);
  Buf<i64> rarea(S), rbase(S), quads(5 * n);
  Buf<i32> nrings(1), rec(2 * (S / 8 + 1));  // <= 2 edges per vertex (4 slots), >= 4 edges per ring
  const MaskView v{bits.p, d_start.p, d_slot.p, d_vrow.p, W.p, H.p, n};
  const Slots sl{succ.p, pred.p, {lab0.p, lab1.p}, {jmp0.p, jmp1.p}, {w0.p, w1.p},
                 rcount.p, rkey.p, rarea.p, rbase.p};
  const int rounds = jump_rounds(max_slots);
  for_each(rows, InitSlots{v, sl, quads.p});
  int in_buf = 0;
  for (int k = 0; k < rounds; ++k, in_buf ^= 1) for_each(rows, LabelRound{v, sl, in_buf});
  const i32* lab = sl.lab[in_buf];
  for_each(rows, Pred{v, sl});
  for_each(rows, RankInit{v, sl, lab});
  in_buf = 0;
  for (int k = 0; k < rounds; ++k, in_buf ^= 1) for_each(rows, RankRound{v, sl, in_buf});
  const i32* w = sl.w[in_buf];
  for_each(rows, RingStats{v, sl, lab, nrings.p, rec.p});
  sync();
  // Order the rings as polygonize does and give each its output offset (host: rings are few).
  const i32 R = nrings.p[0];
  std::vector<i32> order(R);
  std::iota(order.begin(), order.end(), 0);
  auto key = [&](i32 k) { return rkey.p[slot_start[rec.p[2 * k]] + rec.p[2 * k + 1]]; };
  std::sort(order.begin(), order.end(), [&](i32 a, i32 b) {
    return rec.p[2 * a] != rec.p[2 * b] ? rec.p[2 * a] < rec.p[2 * b] : key(a) < key(b);
  });
  const i64 p0 = out.pts_start.back();
  i64 total = 0;
  std::vector<i64> per_mask(n, 0);
  for (i32 k : order) {
    const i64 at = slot_start[rec.p[2 * k]] + rec.p[2 * k + 1];
    rbase.p[at] = total;
    total += rcount.p[at];
    out.pts_start.push_back(p0 + total);
    out.ring_area.push_back(rarea.p[at] / 2);
    per_mask[rec.p[2 * k]] += 1;
  }
  for (i64 k = 0; k < n; ++k) out.ring_start.push_back(out.ring_start.back() + per_mask[k]);
  Buf<i32> pts(2 * total);
  for_each(rows, ScatterCorners{v, sl, lab, w, pts.p});
  sync();
  out.pts.insert(out.pts.end(), pts.p, pts.p + 2 * total);
  out.quads.insert(out.quads.end(), quads.p, quads.p + 5 * n);
}

}  // namespace

// ===================================================================== API
bool cuda_build() {
#ifdef __CUDACC__
  return true;
#else
  return false;
#endif
}

int device_count() {
#ifdef __CUDACC__
  int n = 0;
  return cudaGetDeviceCount(&n) == cudaSuccess ? n : 0;
#else
  return 0;
#endif
}

PolygonResult check_polygons(const Polygons& in) {
  const i64 n = in.n, nr = n ? in.shape_start[n] : 0, nv = nr ? in.ring_start[nr] : 0;
  PolygonResult out;
  out.mask_start.assign(n + 1, 0);
  std::vector<i64> row_start(n + 1, 0);
  for (i64 k = 0; k < n; ++k) {
    out.mask_start[k + 1] = out.mask_start[k] + i64(in.W[k]) * in.H[k];
    row_start[k + 1] = row_start[k] + in.H[k];
  }
  const i64 npx = out.mask_start[n];
  Buf<i32> xy(2 * nv, in.xy), roles(nr, in.roles), W(n, in.W), H(n, in.H);
  Buf<i64> ring_start(nr + 1, in.ring_start), shape_start(n + 1, in.shape_start);
  Buf<i64> d_row(n + 1, row_start.data()), d_mask(n + 1, out.mask_start.data());
  const Shapes s{xy.p, ring_start.p, shape_start.p, d_row.p, d_mask.p, roles.p, W.p, H.p, n, in.G};
  Buf<u8> raster_wide(n), validate_wide(n), eo(npx), nz(npx), fewer(nr);
  Buf<i32> cxy(2 * nv), len(nr), box(4 * nr), counts(nr * kNumDefects), pairs(n * 2 * kPairs);
  Buf<i64> area2(nr), defects(n * kNumDefects);
  const Rings r{cxy.p, len.p, fewer.p, area2.p, box.p, counts.p, pairs.p};
  for_each(n, Precision{s, raster_wide.p, validate_wide.p});
  for_each(row_start[n], RasterRows{s, raster_wide.p, eo.p, nz.p});
  for_each(nr, CleanRings{s, r});
  for_each(nv, EdgePairs{s, r, validate_wide.p});
  for_each(n, ShapeTotals{s, r, validate_wide.p, defects.p});
  sync();
  out.evenodd.assign(eo.p, eo.p + npx);
  out.nonzero.assign(nz.p, nz.p + npx);
  out.defects.assign(defects.p, defects.p + n * kNumDefects);
  return out;
}

MaskResult check_masks(const Masks& in) {
  MaskResult out;
  out.ring_start = {0};
  out.pts_start = {0};
  for (i64 m0 = 0; m0 < in.n;) {  // chunks of at most kChunkSlots slots (a larger mask goes alone)
    i64 m1 = m0, slots = 0;
    while (m1 < in.n) {
      const i64 s = 4 * (i64(in.W[m1]) + 1) * (i64(in.H[m1]) + 1);
      if (m1 > m0 && slots + s > kChunkSlots) break;
      slots += s, ++m1;
    }
    check_mask_chunk(in, m0, m1, out);
    m0 = m1;
  }
  return out;
}

}  // namespace gpucheck
