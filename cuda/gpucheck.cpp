// gpucheck.cpp -- the single-source checker of gpucheck.h. Every algorithm is written once, as __host__ __device__
// functors run by for_each over independent items (rows, rings, edges, shapes); for_each is a plain loop when the
// C++ compiler builds this file and a kernel launch when nvcc does. The other CPU/GPU differences are in the
// "platform" section (memory, copies, streams, four atomic helpers) and the last step of the raster (below).
//
// Polygons (port of cpp/polycore/polycore.hpp, same exact integer arithmetic, int64 or __int128 chosen per shape
// by the same bounds):
//  * raster: one pixel row per item (one warp per row on the GPU). Each edge crossing the row's centre line adds its
//    sign to the winding of the centres left of the crossing; a suffix sum over a 1024-column window gives the
//    winding, so no sort. Only the even-odd mask is made, bit-packed: 32 pixels per uint32 word, LSB = leftmost. On
//    the GPU the warp's lanes share the edges, add their crossings into the window in shared memory, and build each
//    word with __ballot_sync (lane l = column 32 j + l of word j); lane j keeps word j, so a window's 32 words are
//    stored by one coalesced warp store.
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
// Host side (issues CK-1 to CK-3 of experiments/phase2/perf/issues.md): per lane one stream, one device arena and two
// staging buffers (pinned host memory), all made at the first call and grown to the next power of two, never freed
// per call. A polygon call is: pack the inputs into the pinned buffer, one copy in, the kernels on the lane's stream,
// one copy out, one stream synchronize. A large polygon batch is cut into two halves on two lanes, so the second
// half's copy and the host's packing overlap the first half's kernels (C7). One mutex serializes the calls.
//
// References: the headers of cpp/polycore/polycore.hpp and src/objdel/checker/polygonize.py; CUDA C++ Programming
// Guide (https://docs.nvidia.com/cuda/cuda-c-programming-guide/): __int128 in device code (CUDA >= 11.5, "Supported
// Types"), warp vote and shuffle functions (__ballot_sync, __shfl_sync: every lane named in the mask must execute
// them), __syncwarp (memory ordering among the warp's lanes), page-locked host memory and asynchronous concurrent
// execution (cudaMemcpyAsync overlaps kernels only from pinned memory, on a non-default stream).
#include "gpucheck.h"

#include <algorithm>
#include <cstring>
#include <mutex>
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
using polycore::kMaxRings;
using polycore::kNumDefects;
constexpr i32 kPairs = kMaxRings * kMaxRings;

// ===================================================================== platform
#ifdef __CUDACC__
void cuda_check(cudaError_t e) {
  if (e != cudaSuccess) throw std::runtime_error(std::string("CUDA: ") + cudaGetErrorString(e));
}

using Stream = cudaStream_t;

Stream new_stream() {  // non-blocking (no implicit sync with the legacy default stream), lowest priority: the
  cudaStream_t s = nullptr;  // trainer's kernels, queued at the same time, go first
  int lowest = 0, highest = 0;
  cuda_check(cudaDeviceGetStreamPriorityRange(&lowest, &highest));
  cuda_check(cudaStreamCreateWithPriority(&s, cudaStreamNonBlocking, lowest));
  return s;
}

// A failed call (an OOM in Block::get, say) leaves the runtime's last error set, and the next call's
// cudaGetLastError() after a good launch would return it: clear it at the start of every call.
void clear_error() { (void)cudaGetLastError(); }

void* alloc(bool device, size_t bytes) {
  void* p = nullptr;
  cuda_check(device ? cudaMalloc(&p, bytes) : cudaMallocHost(&p, bytes));
  return p;
}

void release(bool device, void* p) { cuda_check(device ? cudaFree(p) : cudaFreeHost(p)); }

void to_device(void* dst, const void* src, size_t bytes, Stream st) {
  if (bytes) cuda_check(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, st));
}

void to_host(void* dst, const void* src, size_t bytes, Stream st) {
  if (bytes) cuda_check(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, st));
}

void zero(void* p, size_t bytes, Stream st) {
  if (bytes) cuda_check(cudaMemsetAsync(p, 0, bytes, st));
}

void wait(Stream st) { cuda_check(cudaStreamSynchronize(st)); }

template <class F>
__global__ void kernel(i64 n, F f) {
  const i64 i = i64(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < n) f(i);
}

template <class F>
void for_each(Stream st, i64 n, const F& f) {
  constexpr int kThreads = 128;
  if (n <= 0) return;
  kernel<<<unsigned((n + kThreads - 1) / kThreads), kThreads, 0, st>>>(n, f);
  cuda_check(cudaGetLastError());
}
#else
struct Stream {};  // the CPU build runs every step at once, in order

Stream new_stream() { return {}; }
void* alloc(bool, size_t bytes) { return ::operator new(bytes); }
void release(bool, void* p) { ::operator delete(p); }
void to_device(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void to_host(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void zero(void* p, size_t bytes, Stream) { if (bytes) std::memset(p, 0, bytes); }
void wait(Stream) {}
void clear_error() {}

template <class F>
void for_each(Stream, i64 n, const F& f) {
  for (i64 i = 0; i < n; ++i) f(i);
}
#endif

// Grow-only buffer: device memory, or host staging memory (pinned on the GPU build). get(bytes) keeps the buffer when
// it is large enough, else replaces it by one of the next power of two up to 64 MiB, exact above (the polygon path
// stays in a few MiB; one large check_masks call must not keep twice its need). Never freed per call.
template <bool kDevice>
struct Block {
  u8* p = nullptr;
  size_t cap = 0;
  u8* get(size_t bytes) {
    if (bytes <= cap) return p;
    size_t c = 4096;
    while (c < bytes && c < (size_t(64) << 20)) c <<= 1;
    if (c < bytes) c = bytes;
    void* old = p;
    p = nullptr, cap = 0;
    if (old) release(kDevice, old);
    p = static_cast<u8*>(alloc(kDevice, c));
    cap = c;
    return p;
  }
};

struct Lane {  // one stream with its device arena and its staging buffers (in: host -> device, out: device -> host)
  Stream stream{};
  Block<true> dev;
  Block<false> in, out;
};

struct Engine {
  std::mutex mu;
  bool ready = false;
  Lane lanes[2];
  void start() {  // called with mu held
    if (ready) return;
    for (Lane& l : lanes) l.stream = new_stream();
    ready = true;
  }
};

// Never destroyed: at process exit the CUDA runtime may be torn down before static objects, so the driver frees all.
Engine& engine() {
  static Engine* e = new Engine;
  return *e;
}

// Waits for both lanes on every exit of a call, so a call that throws leaves no copy in flight from a staging buffer
// the next call refills. Errors are ignored here: the throwing call reports the first one.
struct Drain {
  Engine& e;
  ~Drain() {
#ifdef __CUDACC__
    if (e.ready)
      for (Lane& l : e.lanes) cudaStreamSynchronize(l.stream);
#endif
  }
};

// Arrays of one call placed in one block, each at a 256-byte aligned offset (as cudaMalloc aligns).
struct Layout {
  size_t size = 0;
  template <class T>
  size_t put(i64 n) {
    const size_t at = (size + 255) / 256 * 256;
    size = at + sizeof(T) * size_t(n > 0 ? n : 0);
    return at;
  }
};

template <class T>
T* ptr(u8* base, size_t off) { return reinterpret_cast<T*>(base + off); }

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
  const i64 *ring_start, *shape_start, *row_start, *word_start;
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
constexpr i32 kWin = 1024;  // columns per suffix-sum window = 32 words of 32 pixels

HD inline i32 words_per_row(i32 W) { return (W + 31) / 32; }

// The centre line of one pixel row in polycore's doubled, grid-scaled coordinates.
template <class T>
struct RowLine {
  T sx, sy, g, yc;
  i32 W;
  HD RowLine(const Shapes& s, i64 sh, i32 row)
      : sx(T(2) * s.W[sh]), sy(T(2) * s.H[sh]), g(T(s.G)), yc(T(2 * row + 1) * T(s.G)), W(s.W[sh]) {}

  // Edge a -> b against the line: false when it does not cross (half-open in y); else col = the number of centres
  // strictly left of the crossing, clamped to [0, W] (polycore collect_keys), and sign = +1 for an edge going down.
  HD bool cross(const i32* a, const i32* b, i32& col, i32& sign) const {
    const T ay = T(a[1]) * sy, by = T(b[1]) * sy;
    if (ay == by) return false;  // horizontal edges never count
    const bool down = by > ay;
    const T y0 = down ? ay : by, dy = down ? by - ay : ay - by;
    if (yc < y0 || yc >= y0 + dy) return false;  // half-open: y0 <= yc < y0 + dy
    const T ax = T(a[0]) * sx, bx = T(b[0]) * sx;
    const T x0 = down ? ax : bx, dx = (down ? bx : ax) - x0;
    const T q = floordiv<T>(x0 * dy + (yc - y0) * dx - 1, g * dy);
    const T cnt = floordiv<T>(q - 1, T(2)) + 1;
    col = cnt < 0 ? 0 : (cnt > T(W) ? W : i32(cnt));
    sign = down ? 1 : -1;
    return true;
  }
};

// Crossings of edges first, first + stride, ... of every ring of the shape, for the window [c0, c1): a crossing adds
// its sign to the centres c0 .. col - 1, kept as buf[col - 1 - c0] for the suffix sum. Returns the sum of the signs of
// the crossings right of the whole window (col >= c1).
template <class T>
HD i32 window_crossings(const Shapes& s, i64 sh, const RowLine<T>& L, i32 c0, i32 c1, i32 first, i32 stride,
                        i32* buf) {
  i32 right = 0;
  for (i64 ring = s.shape_start[sh]; ring < s.shape_start[sh + 1]; ++ring) {
    const i64 b = s.ring_start[ring], n = s.ring_start[ring + 1] - b;
    for (i64 k = first; k < n; k += stride) {
      const i64 v0 = b + k, v1 = b + (k + 1 == n ? 0 : k + 1);
      i32 col = 0, sign = 0;
      if (!L.cross(s.xy + 2 * v0, s.xy + 2 * v1, col, sign)) continue;
      if (col >= c1) right += sign;
      else if (col > c0) atomic_add(buf + (col - 1 - c0), sign);
    }
  }
  return right;
}

// Shape, row and first output word of raster row i (rows of all shapes, in order).
HD inline i64 row_of(const Shapes& s, i64 i, i64& sh, i32& row) {
  sh = segment(s.row_start, s.n, i);
  row = i32(i - s.row_start[sh]);
  return s.word_start[sh] + i64(row) * words_per_row(s.W[sh]);
}

#ifdef __CUDACC__
constexpr unsigned kFull = 0xFFFFFFFFu;
constexpr int kRowWarps = 4;  // rows per block of 128 threads; shared memory 4 x 4 KiB

// One warp per row. Every lane runs every warp function below: the row index, the window bounds and the word count
// depend only on the row, so the whole warp takes the same branches around them and kFull names exactly the active
// lanes (a warp whose row is past the end returns as a whole, before any of them).
template <class T>
__device__ void raster_row_warp(const Shapes& s, i64 sh, i32 row, i32* buf, u32* out) {
  const i32 lane = i32(threadIdx.x % 32);
  const RowLine<T> L(s, sh, row);
  for (i32 c0 = 0; c0 < L.W; c0 += kWin) {
    const i32 c1 = mn(L.W, c0 + kWin);
    __syncwarp(kFull);  // the previous window's reads of buf are done
    for (i32 k = lane; k < kWin; k += 32) buf[k] = 0;
    __syncwarp(kFull);
    i32 carry = window_crossings(s, sh, L, c0, c1, lane, 32, buf);
    __syncwarp(kFull);  // every lane's additions to buf are visible to the warp
    for (int off = 16; off > 0; off >>= 1) carry += __shfl_xor_sync(kFull, carry, off);
    const i32 nw = (c1 - c0 + 31) / 32;
    u32 mine = 0;
    for (i32 j = nw - 1; j >= 0; --j) {  // carry = winding contributed by the columns right of word j
      const i32 c = c0 + 32 * j + lane;
      i32 suf = c < c1 ? buf[32 * j + lane] : 0;  // inclusive suffix sum over lanes lane .. 31
      for (int off = 1; off < 32; off <<= 1) {
        const i32 t = __shfl_down_sync(kFull, suf, off);
        if (lane + off < 32) suf += t;
      }
      const u32 word = __ballot_sync(kFull, c < c1 && ((carry + suf) & 1));
      if (lane == j) mine = word;
      carry += __shfl_sync(kFull, suf, 0);
    }
    if (lane < nw) out[c0 / 32 + lane] = mine;  // 32 consecutive words: one coalesced store
  }
}

__global__ void raster_kernel(Shapes s, const u8* wide, u32* words, i64 rows) {
  __shared__ i32 bufs[kRowWarps][kWin];
  const i32 w = i32(threadIdx.x / 32);
  const i64 i = i64(blockIdx.x) * kRowWarps + w;
  if (i >= rows) return;
  i64 sh = 0;
  i32 row = 0;
  u32* out = words + row_of(s, i, sh, row);
  if (wide[sh]) raster_row_warp<i128>(s, sh, row, bufs[w], out);
  else raster_row_warp<i64>(s, sh, row, bufs[w], out);
}

void raster(Stream st, const Shapes& s, const u8* wide, u32* words, i64 rows) {
  if (rows <= 0) return;
  raster_kernel<<<unsigned((rows + kRowWarps - 1) / kRowWarps), 32 * kRowWarps, 0, st>>>(s, wide, words, rows);
  cuda_check(cudaGetLastError());
}
#else
template <class T>
void raster_row(const Shapes& s, i64 sh, i32 row, u32* out) {
  const RowLine<T> L(s, sh, row);
  i32 buf[kWin];
  for (i32 c0 = 0; c0 < L.W; c0 += kWin) {
    const i32 c1 = mn(L.W, c0 + kWin);
    std::fill(buf, buf + (c1 - c0), 0);
    i32 w = window_crossings(s, sh, L, c0, c1, 0, 1, buf);  // winding(c) = sum of the signs right of centre c
    u32 word = 0;
    for (i32 c = c1 - 1; c >= c0; --c) {  // c0 is a multiple of 32: every word is stored at its column 0
      w += buf[c - c0];
      word |= u32(w & 1) << (c & 31);
      if ((c & 31) == 0) out[c / 32] = word, word = 0;
    }
  }
}

void raster(Stream, const Shapes& s, const u8* wide, u32* words, i64 rows) {
  for (i64 i = 0; i < rows; ++i) {
    i64 sh = 0;
    i32 row = 0;
    u32* out = words + row_of(s, i, sh, row);
    if (wide[sh]) raster_row<i128>(s, sh, row, out);
    else raster_row<i64>(s, sh, row, out);
  }
}
#endif

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

struct RingInfo {  // one ring of a mask chunk, gathered for the host's ordering
  i32 m, e, key, count;  // mask (chunk index), start slot, first east edge, corners
  i64 area2;             // twice the signed area
};

struct GatherRings {
  MaskView v;
  Slots sl;
  const i32* rec;
  RingInfo* info;
  HD void operator()(i64 k) const {
    const i32 m = rec[2 * k], e = rec[2 * k + 1];
    const i64 slot = v.slot_start[m] + e;
    info[k] = RingInfo{m, e, sl.rkey[slot], sl.rcount[slot], sl.rarea[slot]};
  }
};

struct SetBase {  // each ring's output offset, at its start slot
  MaskView v;
  Slots sl;
  const i32* rec;
  const i64* base;
  HD void operator()(i64 k) const { sl.rbase[v.slot_start[rec[2 * k]] + rec[2 * k + 1]] = base[k]; }
};

// Masks [m0, m1) on one lane: copy in, rings and corner ranks, then three short round trips to the host (the ring
// count; the rings and bit-quads; the corners after the host ordered the rings).
void check_mask_chunk(const Masks& in, i64 m0, i64 m1, MaskResult& out, Lane& L) {
  const i64 n = m1 - m0;
  std::vector<i64> slot_start(n + 1, 0), vrow_start(n + 1, 0);
  i64 max_slots = 1;
  for (i64 k = 0; k < n; ++k) {
    const i64 W = in.W[m0 + k], H = in.H[m0 + k];
    slot_start[k + 1] = slot_start[k] + 4 * (W + 1) * (H + 1);
    vrow_start[k + 1] = vrow_start[k] + H + 1;
    max_slots = std::max(max_slots, 4 * (W + 1) * (H + 1));
  }
  const i64 npx = in.start[m1] - in.start[m0], S = slot_start[n], rows = vrow_start[n];
  const i64 rmax = S / 8 + 1;  // <= 2 edges per vertex (4 slots), >= 4 edges per ring
  Layout lay;
  const size_t o_bits = lay.put<u8>(npx), o_start = lay.put<i64>(n + 1), o_slot = lay.put<i64>(n + 1),
               o_vrow = lay.put<i64>(n + 1), o_W = lay.put<i32>(n), o_H = lay.put<i32>(n);
  const size_t in_bytes = lay.size;
  size_t o_i32[9];  // succ, pred, lab 0/1, jmp 0/1, w 0/1, rkey
  for (size_t& o : o_i32) o = lay.put<i32>(S);
  const size_t o_rbase = lay.put<i64>(S), o_rec = lay.put<i32>(2 * rmax), o_info = lay.put<RingInfo>(rmax),
               o_base = lay.put<i64>(rmax), o_pts = lay.put<i32>(S);  // <= S / 2 corners
  const size_t o_rcount = lay.put<i32>(S), o_rarea = lay.put<i64>(S), o_quads = lay.put<i64>(5 * n),
               o_nrings = lay.put<i32>(1);  // [o_rcount, end): zeroed
  const size_t info_bytes = (sizeof(RingInfo) * size_t(rmax) + 255) / 256 * 256;
  u8* d = L.dev.get(lay.size);
  u8* h = L.in.get(std::max(in_bytes, sizeof(i64) * size_t(rmax)));
  u8* o = L.out.get(std::max(info_bytes + sizeof(i64) * size_t(5 * n), sizeof(i32) * size_t(S)));
  std::copy_n(in.bits + in.start[m0], npx, ptr<u8>(h, o_bits));
  for (i64 k = 0; k <= n; ++k) ptr<i64>(h, o_start)[k] = in.start[m0 + k] - in.start[m0];
  std::copy_n(slot_start.data(), n + 1, ptr<i64>(h, o_slot));
  std::copy_n(vrow_start.data(), n + 1, ptr<i64>(h, o_vrow));
  std::copy_n(in.W + m0, n, ptr<i32>(h, o_W));
  std::copy_n(in.H + m0, n, ptr<i32>(h, o_H));

  const Stream st = L.stream;
  const MaskView v{ptr<u8>(d, o_bits), ptr<i64>(d, o_start), ptr<i64>(d, o_slot), ptr<i64>(d, o_vrow), ptr<i32>(d, o_W),
                   ptr<i32>(d, o_H), n};
  auto slots = [&](int k) { return ptr<i32>(d, o_i32[k]); };
  const Slots sl{slots(0), slots(1), {slots(2), slots(3)}, {slots(4), slots(5)}, {slots(6), slots(7)},
                 ptr<i32>(d, o_rcount), slots(8), ptr<i64>(d, o_rarea), ptr<i64>(d, o_rbase)};
  i32* rec = ptr<i32>(d, o_rec);
  to_device(d, h, in_bytes, st);
  zero(d + o_rcount, lay.size - o_rcount, st);
  const int rounds = jump_rounds(max_slots);
  for_each(st, rows, InitSlots{v, sl, ptr<i64>(d, o_quads)});
  int in_buf = 0;
  for (int k = 0; k < rounds; ++k, in_buf ^= 1) for_each(st, rows, LabelRound{v, sl, in_buf});
  const i32* lab = sl.lab[in_buf];
  for_each(st, rows, Pred{v, sl});
  for_each(st, rows, RankInit{v, sl, lab});
  in_buf = 0;
  for (int k = 0; k < rounds; ++k, in_buf ^= 1) for_each(st, rows, RankRound{v, sl, in_buf});
  const i32* w = sl.w[in_buf];
  for_each(st, rows, RingStats{v, sl, lab, ptr<i32>(d, o_nrings), rec});
  to_host(o, d + o_nrings, sizeof(i32), st);
  wait(st);
  const i32 R = *ptr<i32>(o, 0);
  if (R < 0 || R > rmax) throw std::logic_error("gpucheck: ring count out of its bound");

  for_each(st, R, GatherRings{v, sl, rec, ptr<RingInfo>(d, o_info)});
  to_host(o, d + o_info, sizeof(RingInfo) * size_t(R), st);
  to_host(o + info_bytes, d + o_quads, sizeof(i64) * size_t(5 * n), st);
  wait(st);
  // Order the rings as polygonize does and give each its output offset (host: rings are few).
  const RingInfo* info = ptr<RingInfo>(o, 0);
  std::vector<i32> order(R);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](i32 a, i32 b) {
    return info[a].m != info[b].m ? info[a].m < info[b].m : info[a].key < info[b].key;
  });
  const i64 p0 = out.pts_start.back();
  i64 total = 0;
  std::vector<i64> per_mask(n, 0);
  i64* base = ptr<i64>(h, 0);  // the copy in is done (wait above): its staging buffer is free
  for (i32 k : order) {
    base[k] = total;
    total += info[k].count;
    out.pts_start.push_back(p0 + total);
    out.ring_area.push_back(info[k].area2 / 2);
    per_mask[info[k].m] += 1;
  }
  for (i64 k = 0; k < n; ++k) out.ring_start.push_back(out.ring_start.back() + per_mask[k]);
  const i64* quads = ptr<i64>(o + info_bytes, 0);
  out.quads.insert(out.quads.end(), quads, quads + 5 * n);
  if (2 * total > S) throw std::logic_error("gpucheck: corner count out of its bound");

  to_device(d + o_base, base, sizeof(i64) * size_t(R), st);
  for_each(st, R, SetBase{v, sl, rec, ptr<i64>(d, o_base)});
  for_each(st, rows, ScatterCorners{v, sl, lab, w, ptr<i32>(d, o_pts)});
  to_host(o, d + o_pts, sizeof(i32) * size_t(2 * total), st);
  wait(st);
  out.pts.insert(out.pts.end(), ptr<i32>(o, 0), ptr<i32>(o, 0) + 2 * total);
}

// ---------------------------------------------------------------- polygon calls
constexpr i64 kSplitMin = 64;  // polygon batches of at least this many shapes run as two halves on the two lanes

// One half (or all) of a polygon call: shapes [a, b), their arrays at offsets of one lane's arena. The staging buffers
// mirror it: [0, in_bytes) is the copy in, [defects, size) the copy out.
struct PolyChunk {
  i64 a, b, n, r0, nr, v0, nv, rows, nwords;
  size_t xy, roles, W, H, ring_start, shape_start, row_start, word_start, in_bytes;
  size_t raster_wide, validate_wide, fewer, cxy, len, box, area2, counts, pairs, defects, words, size;
};

PolyChunk plan(const Polygons& in, const std::vector<i64>& word_start, i64 a, i64 b) {
  PolyChunk c{};
  c.a = a, c.b = b, c.n = b - a;
  c.r0 = in.shape_start[a], c.nr = in.shape_start[b] - c.r0;
  c.v0 = in.ring_start[c.r0], c.nv = in.ring_start[c.r0 + c.nr] - c.v0;
  c.rows = 0;
  for (i64 s = a; s < b; ++s) c.rows += in.H[s];
  c.nwords = word_start[b] - word_start[a];
  Layout l;
  c.xy = l.put<i32>(2 * c.nv), c.roles = l.put<i32>(c.nr), c.W = l.put<i32>(c.n), c.H = l.put<i32>(c.n);
  c.ring_start = l.put<i64>(c.nr + 1), c.shape_start = l.put<i64>(c.n + 1), c.row_start = l.put<i64>(c.n + 1);
  c.word_start = l.put<i64>(c.n + 1);
  c.in_bytes = l.size;
  c.raster_wide = l.put<u8>(c.n), c.validate_wide = l.put<u8>(c.n), c.fewer = l.put<u8>(c.nr);
  c.cxy = l.put<i32>(2 * c.nv), c.len = l.put<i32>(c.nr), c.box = l.put<i32>(4 * c.nr), c.area2 = l.put<i64>(c.nr);
  c.counts = l.put<i32>(c.nr * kNumDefects), c.pairs = l.put<i32>(c.n * 2 * kPairs);  // zeroed, with defects
  c.defects = l.put<i64>(c.n * kNumDefects), c.words = l.put<u32>(c.nwords);
  c.size = l.size;
  return c;
}

void pack(const Polygons& in, const std::vector<i64>& word_start, const PolyChunk& c, u8* h) {
  std::copy_n(in.xy + 2 * c.v0, 2 * c.nv, ptr<i32>(h, c.xy));  // copy_n: a count of 0 never reads (null data)
  std::copy_n(in.roles + c.r0, c.nr, ptr<i32>(h, c.roles));
  std::copy_n(in.W + c.a, c.n, ptr<i32>(h, c.W));
  std::copy_n(in.H + c.a, c.n, ptr<i32>(h, c.H));
  for (i64 r = 0; r <= c.nr; ++r) ptr<i64>(h, c.ring_start)[r] = in.ring_start[c.r0 + r] - c.v0;
  i64* row_start = ptr<i64>(h, c.row_start);
  row_start[0] = 0;
  for (i64 s = 0; s < c.n; ++s) {
    ptr<i64>(h, c.shape_start)[s] = in.shape_start[c.a + s] - c.r0;
    ptr<i64>(h, c.word_start)[s] = word_start[c.a + s] - word_start[c.a];
    row_start[s + 1] = row_start[s] + in.H[c.a + s];
  }
  ptr<i64>(h, c.shape_start)[c.n] = c.nr;
  ptr<i64>(h, c.word_start)[c.n] = c.nwords;
}

// Queues the chunk on its lane: copy in, zero the counters, the kernels, copy out. Nothing waits.
void launch(const PolyChunk& c, i64 G, Lane& L) {
  u8* d = L.dev.p;
  const Stream st = L.stream;
  const Shapes s{ptr<i32>(d, c.xy), ptr<i64>(d, c.ring_start), ptr<i64>(d, c.shape_start), ptr<i64>(d, c.row_start),
                 ptr<i64>(d, c.word_start), ptr<i32>(d, c.roles), ptr<i32>(d, c.W), ptr<i32>(d, c.H), c.n, G};
  const Rings r{ptr<i32>(d, c.cxy), ptr<i32>(d, c.len), ptr<u8>(d, c.fewer), ptr<i64>(d, c.area2), ptr<i32>(d, c.box),
                ptr<i32>(d, c.counts), ptr<i32>(d, c.pairs)};
  u8 *raster_wide = ptr<u8>(d, c.raster_wide), *validate_wide = ptr<u8>(d, c.validate_wide);
  to_device(d, L.in.p, c.in_bytes, st);
  zero(d + c.counts, c.words - c.counts, st);  // counts, pairs, defects: the masks are written whole
  for_each(st, c.n, Precision{s, raster_wide, validate_wide});
  raster(st, s, raster_wide, ptr<u32>(d, c.words), c.rows);
  for_each(st, c.nr, CleanRings{s, r});
  for_each(st, c.nv, EdgePairs{s, r, validate_wide});
  for_each(st, c.n, ShapeTotals{s, r, validate_wide, ptr<i64>(d, c.defects)});
  to_host(L.out.p, d + c.defects, c.size - c.defects, st);
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
  const i64 n = in.n;
  PolygonResult out;
  out.word_start.assign(n + 1, 0);
  for (i64 k = 0; k < n; ++k) out.word_start[k + 1] = out.word_start[k] + i64(in.H[k]) * words_per_row(in.W[k]);
  if (n == 0) return out;
  out.evenodd.resize(size_t(out.word_start[n]));
  out.defects.resize(size_t(n * kNumDefects));
  const i64 half = n >= kSplitMin ? n / 2 : n;
  const PolyChunk chunks[2] = {plan(in, out.word_start, 0, half), plan(in, out.word_start, half, n)};
  const int used = half < n ? 2 : 1;
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  e.start();
  for (int k = 0; k < used; ++k) {  // every buffer is grown before any work is queued
    Lane& L = e.lanes[k];
    L.dev.get(chunks[k].size), L.in.get(chunks[k].in_bytes), L.out.get(chunks[k].size - chunks[k].defects);
  }
  for (int k = 0; k < used; ++k) {  // lane 1 is packed while lane 0 runs
    pack(in, out.word_start, chunks[k], e.lanes[k].in.p);
    launch(chunks[k], in.G, e.lanes[k]);
  }
  for (int k = 0; k < used; ++k) {
    const PolyChunk& c = chunks[k];
    Lane& L = e.lanes[k];
    wait(L.stream);
    std::copy_n(ptr<i64>(L.out.p, 0), c.n * kNumDefects, out.defects.data() + c.a * kNumDefects);
    std::copy_n(ptr<u32>(L.out.p, c.words - c.defects), c.nwords, out.evenodd.data() + out.word_start[c.a]);
  }
  return out;
}

MaskResult check_masks(const Masks& in) {
  MaskResult out;
  out.ring_start = {0};
  out.pts_start = {0};
  if (in.n == 0) return out;
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  e.start();
  for (i64 m0 = 0; m0 < in.n;) {  // chunks of at most kChunkSlots slots (a larger mask goes alone)
    i64 m1 = m0, slots = 0;
    while (m1 < in.n) {
      const i64 s = 4 * (i64(in.W[m1]) + 1) * (i64(in.H[m1]) + 1);
      if (m1 > m0 && slots + s > kChunkSlots) break;
      slots += s, ++m1;
    }
    check_mask_chunk(in, m0, m1, out, e.lanes[0]);
    m0 = m1;
  }
  return out;
}

}  // namespace gpucheck
