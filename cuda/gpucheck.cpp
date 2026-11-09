// gpucheck.cpp -- the single-source checker of gpucheck.h. Every algorithm is written once, as __host__ __device__
// functors run by for_each over independent items (rows, rings, edges, shapes); for_each is a plain loop when the
// C++ compiler builds this file and a kernel launch when nvcc does. The other CPU/GPU differences are in the
// "platform" section (memory, copies, streams, four atomic helpers) and the last step of the raster (below).
//
// Polygons (port of the archived CPU checker, archive/cpu_checker/cpp/polycore/polycore.hpp: same exact integer
// arithmetic, int64 or __int128 chosen per shape by the same bounds):
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
// References: the headers of archive/cpu_checker/cpp/polycore/polycore.hpp and src/objdel/checker/polygonize.py; CUDA C++ Programming
// Guide (https://docs.nvidia.com/cuda/cuda-c-programming-guide/): __int128 in device code (CUDA >= 11.5, "Supported
// Types"), warp vote and shuffle functions (__ballot_sync, __shfl_sync: every lane named in the mask must execute
// them), __syncwarp (memory ordering among the warp's lanes), page-locked host memory and asynchronous concurrent
// execution (cudaMemcpyAsync overlaps kernels only from pinned memory, on a non-default stream).
#include "gpucheck.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#include <nvtx3/nvtx3.hpp>
#define HD __host__ __device__
#else
#define HD
#endif

// NVTX ranges for nsys, in the "gpucheck" domain (CUDA build only; they cost a few instructions without a profiler).
// Each use declares its own variable, so two ranges can share a scope; a range ends with its scope.
#define GC_CAT2(a, b) a##b
#define GC_CAT(a, b) GC_CAT2(a, b)
#ifdef __CUDACC__
struct gc_domain {
  static constexpr char const* name{"gpucheck"};
};
#define GC_RANGE(name) nvtx3::scoped_range_in<gc_domain> GC_CAT(gc_range_, __LINE__) { name }
#else
#define GC_RANGE(name)
#endif

namespace gpucheck {
namespace {

using i128 = __int128;
constexpr i32 kPairs = kMaxRings * kMaxRings;

// The engine counters (gpucheck.h Stats). Written with the engine mutex held, read under it (stats()).
Stats& counters() {
  static Stats s;
  return s;
}

// ===================================================================== platform
// Every host wait, stream sync, allocation and release goes through the helpers of this section, which count them.
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

void* alloc(bool device, size_t bytes) {  // cudaMalloc and cudaMallocHost synchronize the device
  void* p = nullptr;
  cuda_check(device ? cudaMalloc(&p, bytes) : cudaMallocHost(&p, bytes));
  ++counters().sync_alloc;
  if (device) counters().live_device_bytes += i64(bytes);
  return p;
}

void release(bool device, void* p, size_t bytes) {
  ++counters().sync_free;
  if (device) counters().live_device_bytes -= i64(bytes);
  cuda_check(device ? cudaFree(p) : cudaFreeHost(p));
}

void to_device(void* dst, const void* src, size_t bytes, Stream st) {
  if (bytes) cuda_check(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, st));
}

void to_host(void* dst, const void* src, size_t bytes, Stream st) {
  if (bytes) cuda_check(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, st));
}

void zero(void* p, size_t bytes, Stream st) {
  if (bytes) cuda_check(cudaMemsetAsync(p, 0, bytes, st));
}

void wait(Stream st) {
  ++counters().waits;
  cuda_check(cudaStreamSynchronize(st));
}

void drain_sync(Stream st) {  // errors ignored: the call that throws reports the first one
  ++counters().drain_syncs;
  (void)cudaStreamSynchronize(st);
}

// Stage marks (Timer): CUDA events on the lane's stream.
using Mark = cudaEvent_t;
Mark mark(Stream st) {
  cudaEvent_t ev = nullptr;
  cuda_check(cudaEventCreate(&ev));
  cuda_check(cudaEventRecord(ev, st));
  return ev;
}
double elapsed_ms(Mark a, Mark b) {
  float ms = 0;
  cuda_check(cudaEventElapsedTime(&ms, a, b));
  return ms;
}
void unmark(Mark m) { (void)cudaEventDestroy(m); }

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
// Counted as the CUDA build counts them after the stream-ordered pool (D5c): device memory as pool allocations, host
// staging as synchronizing ones.
void* alloc(bool device, size_t bytes) {
  void* p = ::operator new(bytes);
  ++(device ? counters().pool_allocs : counters().sync_alloc);
  if (device) counters().live_device_bytes += i64(bytes);
  return p;
}
void release(bool device, void* p, size_t bytes) {
  ++(device ? counters().pool_frees : counters().sync_free);
  if (device) counters().live_device_bytes -= i64(bytes);
  ::operator delete(p);
}
void to_device(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void to_host(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void zero(void* p, size_t bytes, Stream) { if (bytes) std::memset(p, 0, bytes); }
void wait(Stream) { ++counters().waits; }
void drain_sync(Stream) { ++counters().drain_syncs; }
void clear_error() {}

// Stage marks (Timer): the host clock, since every step has run when the next one is queued.
using Mark = std::chrono::steady_clock::time_point;
Mark mark(Stream) { return std::chrono::steady_clock::now(); }
double elapsed_ms(Mark a, Mark b) { return std::chrono::duration<double, std::milli>(b - a).count(); }
void unmark(Mark) {}

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
    const size_t old_cap = cap;
    p = nullptr, cap = 0;
    if (kDevice) counters().arena_bytes -= i64(old_cap);
    if (old) release(kDevice, old, old_cap);
    p = static_cast<u8*>(alloc(kDevice, c));
    cap = c;
    if (kDevice) counters().arena_bytes += i64(cap);
    return p;
  }
};

struct Lane {  // one stream with its device arenas (dev: checker, score: scorer) and its staging buffers (in: host ->
  Stream stream{};  // device, out: device -> host)
  Block<true> dev, score;
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
#ifdef __CUDACC__
    counters().contexts = 1;  // the streams live in the device's primary context
#endif
  }
};

// Never destroyed: at process exit the CUDA runtime may be torn down before static objects, so the driver frees all.
Engine& engine() {
  static Engine* e = new Engine;
  return *e;
}

// Waits for both lanes on every exit of a call, so a call that throws leaves no copy in flight from a staging buffer
// the next call refills.
struct Drain {
  Engine& e;
  ~Drain() {
    if (e.ready)
      for (Lane& l : e.lanes) drain_sync(l.stream);
  }
};

// The per-stage times of one call (gpucheck.h Times): a mark before and after each stage on its lane's stream, read
// by read() after the call's last wait, which covers every mark. Without `out`, a stage is only its NVTX range.
struct Timer {
  struct Span {
    const char* name;
    Mark a, b;
  };
  Times* out;
  std::vector<Span> spans;
  template <class F>
  void stage(const char* name, Stream st, const F& f) {
    GC_RANGE(name);
    if (!out) return f();
    const Mark a = mark(st);
    f();
    spans.push_back(Span{name, a, mark(st)});
  }
  void read() {
    for (const Span& s : spans) {
      auto it = std::find_if(out->begin(), out->end(), [&](const auto& t) { return t.first == s.name; });
      if (it == out->end()) it = out->insert(out->end(), {s.name, 0.0});
      it->second += elapsed_ms(s.a, s.b);
    }
  }
  ~Timer() {
    for (const Span& s : spans) unmark(s.a), unmark(s.b);
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

HD inline i32 atomic_min_old(i32* p, i32 v) {  // returns the old value
#ifdef __CUDA_ARCH__
  return atomicMin(p, v);
#else
  const i32 old = *p;
  if (v < old) *p = v;
  return old;
#endif
}

HD inline i32 load(const i32* p) {  // a fresh read of memory other threads write (union-find parents)
#ifdef __CUDA_ARCH__
  return *static_cast<const volatile i32*>(p);
#else
  return *p;
#endif
}

HD inline int popc(u32 x) {
#ifdef __CUDA_ARCH__
  return __popc(x);
#else
  return __builtin_popcount(x);
#endif
}

HD inline int ctz(u32 x) {  // x != 0
#ifdef __CUDA_ARCH__
  return __ffs(x) - 1;
#else
  return __builtin_ctz(x);
#endif
}

HD inline int clz(u32 x) {  // x != 0
#ifdef __CUDA_ARCH__
  return __clz(x);
#else
  return __builtin_clz(x);
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
    cnt[kOutOfRange] = oor;
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
    if (distinct < 3) cnt[kFewerThan3Points] = 1;
    else if (a2 == 0) cnt[kZeroArea] = 1;
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
          const int d = kind == kCross ? kSelfCrossing
                      : kind == kPinch ? kVertexPinch
                      : kind == kOverlap ? kCollinearOverlap : kSpike;
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
        at(a, kRingRingCrossing) += pcross[a * kMaxRings + b];
        at(a, kRingRingTouching) += ptouch[a * kMaxRings + b];
      }
    for (i32 a = 0; a < R; ++a) {
      if (fewer(a)) continue;
      const i64 a2 = r.area2[r0 + a];
      if ((roles[a] < 0 && a2 < 0) || (roles[a] >= 0 && a2 > 0)) at(a, kWrongWinding) = 1;
    }
    for (i32 h = 0; h < R; ++h) {
      const i32 p = roles[h];
      if (p < 0 || fewer(h) || fewer(p)) continue;
      if (any_witness<T>(r0, R, h, -1, apart(h, p), p, false)) at(h, kHoleOutsideOuter) = 1;
    }
    for (i32 a = 0; a < R; ++a) {
      if (fewer(a)) continue;
      for (i32 b = a + 1; b < R; ++b) {
        if (fewer(b)) continue;
        if (roles[a] >= 0 && roles[a] == roles[b]) {  // sibling holes
          const bool single = apart(a, b);
          if (any_witness<T>(r0, R, b, 1, single, a, false) || any_witness<T>(r0, R, a, 1, single, b, false)) {
            at(a, kNestedHoles) += 1, at(b, kNestedHoles) += 1;
            tot[kNestedHoles] += 1;
          }
        } else if (roles[a] < 0 && roles[b] < 0) {  // two outers
          if (any_witness<T>(r0, R, b, 1, apart_region(b, a), a, true) ||
              any_witness<T>(r0, R, a, 1, apart_region(a, b), b, true)) {
            at(a, kNestedOuters) += 1, at(b, kNestedOuters) += 1;
            tot[kNestedOuters] += 1;
          }
        }
      }
    }
    for (i32 a = 0; a < R; ++a)
      for (int d = 0; d < kNumDefects; ++d)
        if (d != kRingRingCrossing && d != kRingRingTouching && d != kNestedHoles &&
            d != kNestedOuters)
          tot[d] += at(a, d);
    for (i32 a = 0; a < R; ++a)
      for (i32 b = a + 1; b < R; ++b) {
        tot[kRingRingCrossing] += pcross[a * kMaxRings + b];
        tot[kRingRingTouching] += ptouch[a * kMaxRings + b];
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
void check_mask_chunk(const Masks& in, i64 m0, i64 m1, MaskResult& out, Lane& L, Timer& t) {
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
  const i32 *lab = nullptr, *w = nullptr;
  t.stage("ccl", st, [&] {  // ring labels: pointer jumping to each cycle's minimum slot
    for_each(st, rows, InitSlots{v, sl, ptr<i64>(d, o_quads)});
    int in_buf = 0;
    for (int k = 0; k < rounds; ++k, in_buf ^= 1) for_each(st, rows, LabelRound{v, sl, in_buf});
    lab = sl.lab[in_buf];
  });
  t.stage("rings", st, [&] {  // corner ranks and the ring records
    for_each(st, rows, Pred{v, sl});
    for_each(st, rows, RankInit{v, sl, lab});
    int in_buf = 0;
    for (int k = 0; k < rounds; ++k, in_buf ^= 1) for_each(st, rows, RankRound{v, sl, in_buf});
    w = sl.w[in_buf];
    for_each(st, rows, RingStats{v, sl, lab, ptr<i32>(d, o_nrings), rec});
  });
  to_host(o, d + o_nrings, sizeof(i32), st);
  wait(st);
  const i32 R = *ptr<i32>(o, 0);
  if (R < 0 || R > rmax) throw std::logic_error("gpucheck: ring count out of its bound");

  t.stage("rings", st, [&] { for_each(st, R, GatherRings{v, sl, rec, ptr<RingInfo>(d, o_info)}); });
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
  t.stage("trace", st, [&] {  // every corner to its output position
    for_each(st, R, SetBase{v, sl, rec, ptr<i64>(d, o_base)});
    for_each(st, rows, ScatterCorners{v, sl, lab, w, ptr<i32>(d, o_pts)});
  });
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

// Queues the chunk on its lane: copy in, zero the counters, the kernels, copy out (the defects, and the masks when
// `masks`; else they stay in the lane's arena for the scorer). Nothing waits.
void launch(const PolyChunk& c, i64 G, Lane& L, bool masks, Timer& t) {
  u8* d = L.dev.p;
  const Stream st = L.stream;
  const Shapes s{ptr<i32>(d, c.xy), ptr<i64>(d, c.ring_start), ptr<i64>(d, c.shape_start), ptr<i64>(d, c.row_start),
                 ptr<i64>(d, c.word_start), ptr<i32>(d, c.roles), ptr<i32>(d, c.W), ptr<i32>(d, c.H), c.n, G};
  const Rings r{ptr<i32>(d, c.cxy), ptr<i32>(d, c.len), ptr<u8>(d, c.fewer), ptr<i64>(d, c.area2), ptr<i32>(d, c.box),
                ptr<i32>(d, c.counts), ptr<i32>(d, c.pairs)};
  u8 *raster_wide = ptr<u8>(d, c.raster_wide), *validate_wide = ptr<u8>(d, c.validate_wide);
  to_device(d, L.in.p, c.in_bytes, st);
  zero(d + c.counts, c.words - c.counts, st);  // counts, pairs, defects: the masks are written whole
  t.stage("raster", st, [&] {
    for_each(st, c.n, Precision{s, raster_wide, validate_wide});
    raster(st, s, raster_wide, ptr<u32>(d, c.words), c.rows);
  });
  t.stage("validity", st, [&] {
    for_each(st, c.nr, CleanRings{s, r});
    for_each(st, c.nv, EdgePairs{s, r, validate_wide});
    for_each(st, c.n, ShapeTotals{s, r, validate_wide, ptr<i64>(d, c.defects)});
  });
  to_host(L.out.p, d + c.defects, (masks ? c.size : c.words) - c.defects, st);
}

// ===================================================================== scoring
// The reward's per-object terms (objdel.reward.score), computed where the mask is: the masks of a call stay in the
// lane's arena, every expert table is resident (add_tables), and only per-object numbers are copied back. Bit-identical
// to the host scorer it replaced (docs/research/notes/gpu_scoring.md):
//  * counts (area, intersections with each rater and the vote) are integers; IoU = inter / union in float64;
//  * band depth D sums float32 distances that are 0 or >= 1, hence multiples of 2^-23: each sum is an exact integer
//    count of 2^-23 (per row in int64, rows in order in int128), so its float64 value equals that of any float64
//    summation order (numpy's) while the sum is below 2^30; F5, g_area and the shares are the same IEEE operations in
//    the same order as the Python code (built without FMA contraction, CMakeLists.txt);
//  * the signature: union-find labeling on the mask's bounding box padded by one pixel (8-connected foreground,
//    4-connected background; the outside is the background component of pixel 0), the steps of
//    objdel.checker.topology.signature: drop components below min_area, fill holes below min_area, relabel, count each
//    hole for the component above its first raster pixel, sort descending.
// One call: the box of every mask (round trip 1: the crop layout), the labeling, row partials and per-object terms
// (round trip 2: the terms and signature lengths), the signatures (round trip 3).

struct TableView {  // device view of resident tables: bit planes of H rows of words_per_row(W) words, LSB = leftmost
  const u32 *raters, *uni, *core, *vote;  // raters: k planes
  const float* dist;
  const i64 *rater_area, *sig_start;
  const i32* sig;
  const u8* voter;
  i64 n_vote, perimeter;
  double t5;
  i32 H, W, k, min_area;
};

struct Registry {  // resident tables; handle = index (reused once dropped)
  std::vector<u8*> block;
  std::vector<size_t> bytes;
  std::vector<TableView> view;
  std::vector<i64> free;
  i64 live = 0;
};

Registry& registry() {  // never destroyed, as engine()
  static Registry* r = new Registry;
  return *r;
}

void pack_plane(const u8* px, i32 H, i32 W, u32* out) {  // out zeroed
  const i64 wpr = words_per_row(W);
  for (i32 y = 0; y < H; ++y)
    for (i32 x = 0; x < W; ++x)
      if (px[i64(y) * W + x]) out[y * wpr + x / 32] |= u32(1) << (x % 32);
}

std::vector<TableView> views(const i64* handles, i64 n, const i32* W, const i32* H) {  // called with the mutex held
  const Registry& r = registry();
  std::vector<TableView> v(static_cast<size_t>(n));
  for (i64 s = 0; s < n; ++s) {
    const i64 h = handles[s];
    if (h < 0 || h >= i64(r.block.size()) || !r.block[size_t(h)]) throw std::invalid_argument("unknown tables handle");
    v[size_t(s)] = r.view[size_t(h)];
    if (v[size_t(s)].W != W[s] || v[size_t(s)].H != H[s])
      throw std::invalid_argument("an object's frame must be its tables' frame");
  }
  return v;
}

constexpr int kPart = 5;  // row partials after the k rater intersections: vote intersection, spill count and sum, miss
                          // count and sum (sums in units of 2^-23)

struct Obj {  // one object of a scoring call
  TableView t;
  i64 words;           // its mask: t.H rows of words_per_row(t.W) words
  i64 part;            // its row partials: t.H rows of t.k + kPart
  i64 crop;            // its crop pixels: bh + 2 rows of bw + 2 (none for an empty mask)
  i64 area;
  i32 bx, by, bw, bh;  // the mask's bounding box
};

struct ObjOut {  // copied back per object
  i64 area, ncomp;
  double iou_tc, f5, g_area, topo_freq;
};

// Per frame row: first and last foreground column, foreground pixels.
struct RowBox {
  const Obj* obj;
  const i64* row_start;
  i64 n;
  const u32* words;
  i32* rbox;
  HD void operator()(i64 i) const {
    const i64 o = segment(row_start, n, i);
    const i32 wpr = words_per_row(obj[o].t.W);
    const u32* w = words + obj[o].words + (i - row_start[o]) * wpr;
    i32 first = -1, last = -1, count = 0;
    for (i32 j = 0; j < wpr; ++j) {
      if (!w[j]) continue;
      if (first < 0) first = 32 * j + ctz(w[j]);
      last = 32 * j + 31 - clz(w[j]);
      count += popc(w[j]);
    }
    rbox[3 * i] = first, rbox[3 * i + 1] = last, rbox[3 * i + 2] = count;
  }
};

// Per object: box (x, y, w, h) and area.
struct ObjBox {
  const i64* row_start;
  const i32* rbox;
  i64* box;
  HD void operator()(i64 o) const {
    i64 area = 0;
    i32 x0 = 0, x1 = 0, y0 = -1, y1 = 0;
    for (i64 i = row_start[o]; i < row_start[o + 1]; ++i) {
      const i32* r = rbox + 3 * i;
      if (!r[2]) continue;
      x0 = y0 < 0 ? r[0] : mn(x0, r[0]), x1 = y0 < 0 ? r[1] : mx(x1, r[1]);
      y1 = i32(i - row_start[o]), area += r[2];
      if (y0 < 0) y0 = y1;
    }
    i64* b = box + 5 * o;
    b[0] = area ? x0 : 0, b[1] = area ? y0 : 0, b[2] = area ? x1 - x0 + 1 : 0, b[3] = area ? y1 - y0 + 1 : 0;
    b[4] = area;
  }
};

// ---------------------------------------------------------------- signature: union-find labeling on the crops
struct Crops {  // the call's crops, one item per crop row; pixel p of a crop = y (bw + 2) + x, local to it
  const Obj* obj;
  const i64* crow_start;
  i64 n;
  u8* img;             // 1 = foreground
  i32 *fg, *bg, *aux;  // labels (root = the component's first raster pixel), per-root counts
  i32* ncomp;          // per object: components after the filter
  HD i64 row(i64 r, i32& y) const {
    const i64 o = segment(crow_start, n, r);
    y = i32(r - crow_start[o]);
    return o;
  }
};

HD inline i32 find(const i32* par, i32 x) {
  for (i32 p = load(par + x); p != x; p = load(par + x)) x = p;
  return x;
}

// Union by linking the larger root under the smaller (Playne and Hawick 2018: a failed atomicMin means the root was
// linked meanwhile; retry from the value it had). Roots only ever get smaller parents, so the final root of a
// component is its smallest pixel index, its first raster pixel.
HD inline void unite(i32* par, i32 a, i32 b) {
  for (;;) {
    a = find(par, a), b = find(par, b);
    if (a == b) return;
    if (a < b) {
      const i32 t = a;
      a = b, b = t;
    }
    const i32 old = atomic_min_old(par + a, b);
    if (old == a) return;
    a = old;
  }
}

struct CropFill {
  Crops c;
  const u32* words;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Obj& b = c.obj[c.row(r, y)];
    const i32 cw = b.bw + 2, wpr = words_per_row(b.t.W);
    u8* img = c.img + b.crop + i64(y) * cw;
    const bool inside = y >= 1 && y <= b.bh;
    const u32* row = words + b.words + i64(b.by + (inside ? y - 1 : 0)) * wpr;
    for (i32 x = 0; x < cw; ++x) {
      const i32 fx = b.bx + x - 1;
      img[x] = (inside && x >= 1 && x <= b.bw) ? u8(row[fx / 32] >> (fx % 32) & 1) : u8(0);
    }
  }
};

// The three labeling steps of one class: which = 1 (foreground, 8-connected) or 0 (background, 4-connected).
struct LabelInit {
  Crops c;
  int which;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Obj& b = c.obj[c.row(r, y)];
    const i32 cw = b.bw + 2;
    i32* par = (which ? c.fg : c.bg) + b.crop;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p) par[p] = c.img[b.crop + p] == which ? p : -1;
  }
};

struct LabelMerge {
  Crops c;
  int which;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Obj& b = c.obj[c.row(r, y)];
    const i32 cw = b.bw + 2;
    const u8* img = c.img + b.crop;
    i32* par = (which ? c.fg : c.bg) + b.crop;
    for (i32 x = 0; x < cw; ++x) {
      const i32 p = y * cw + x;
      if (img[p] != which) continue;
      if (x > 0 && img[p - 1] == which) unite(par, p, p - 1);
      if (y == 0) continue;
      if (img[p - cw] == which) unite(par, p, p - cw);
      if (!which) continue;
      if (x > 0 && img[p - cw - 1]) unite(par, p, p - cw - 1);
      if (x + 1 < cw && img[p - cw + 1]) unite(par, p, p - cw + 1);
    }
  }
};

struct LabelFlatten {
  Crops c;
  int which;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Obj& b = c.obj[c.row(r, y)];
    const i32 cw = b.bw + 2;
    i32* par = (which ? c.fg : c.bg) + b.crop;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p)
      if (c.img[b.crop + p] == which) par[p] = find(par, p);
  }
};

struct Count {  // aux[root] = pixels of the class's component (aux zeroed)
  Crops c;
  int which;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Obj& b = c.obj[c.row(r, y)];
    const i32 cw = b.bw + 2;
    const i32* par = (which ? c.fg : c.bg) + b.crop;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p)
      if (c.img[b.crop + p] == which) atomic_add(c.aux + b.crop + par[p], 1);
  }
};

// Size filter (topology.size_filter): which = 1 drops foreground components below min_area, which = 0 fills holes
// (background components other than the outside, root 0) below min_area. aux holds the class's counts.
struct Filter {
  Crops c;
  int which;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Obj& b = c.obj[c.row(r, y)];
    const i32 cw = b.bw + 2;
    const i32* par = (which ? c.fg : c.bg) + b.crop;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p) {
      if (c.img[b.crop + p] != which || (!which && par[p] == 0)) continue;
      if (c.aux[b.crop + par[p]] < b.t.min_area) c.img[b.crop + p] = u8(1 - which);
    }
  }
};

// A hole's first raster pixel is its root; the pixel above it is foreground of the enclosing component (an island
// inside the hole lies strictly below the hole's top row): aux[foreground root] += 1 (aux zeroed).
struct Holes {
  Crops c;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Obj& b = c.obj[c.row(r, y)];
    const i32 cw = b.bw + 2;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p)
      if (p != 0 && c.img[b.crop + p] == 0 && c.bg[b.crop + p] == p)
        atomic_add(c.aux + b.crop + c.fg[b.crop + p - cw], 1);
  }
};

// Every component's hole count into the object's list: the first ncomp entries of its crop in bg (no longer needed).
struct Gather {
  Crops c;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const i64 o = c.row(r, y);
    const Obj& b = c.obj[o];
    const i32 cw = b.bw + 2;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p)
      if (c.img[b.crop + p] == 1 && c.fg[b.crop + p] == p) c.bg[b.crop + atomic_inc(c.ncomp + o)] = c.aux[b.crop + p];
  }
};

// ---------------------------------------------------------------- terms
HD inline i64 units(float d) { return i64(double(d) * 8388608.0); }  // d in units of 2^-23: exact for 0 and d >= 1

// Per frame row: intersections with each rater and the vote, spill (mask outside the union) and core miss (core outside
// the mask) pixel counts and distance sums. The bits past W are 0 in every plane and mask.
struct ScoreRows {
  const Obj* obj;
  const i64* row_start;
  i64 n;
  const u32* words;
  i64* parts;
  HD void operator()(i64 i) const {
    const i64 o = segment(row_start, n, i);
    const Obj& b = obj[o];
    const TableView& t = b.t;
    const i64 y = i - row_start[o], wpr = words_per_row(t.W);
    const u32* m = words + b.words + y * wpr;
    i64* out = parts + b.part + y * (t.k + kPart);
    for (i32 j = 0; j < t.k; ++j) {
      const u32* r = t.raters + (j * i64(t.H) + y) * wpr;
      i64 s = 0;
      for (i64 w = 0; w < wpr; ++w) s += popc(m[w] & r[w]);
      out[j] = s;
    }
    i64 vote = 0, spill_n = 0, spill = 0, miss_n = 0, miss = 0;
    const float* d = t.dist + y * t.W;
    for (i64 w = 0; w < wpr; ++w) {
      vote += popc(m[w] & t.vote[y * wpr + w]);
      u32 s = m[w] & ~t.uni[y * wpr + w], c = t.core[y * wpr + w] & ~m[w];
      spill_n += popc(s), miss_n += popc(c);
      for (; s; s &= s - 1) spill += units(d[32 * w + ctz(s)]);
      for (; c; c &= c - 1) miss += units(d[32 * w + ctz(c)]);
    }
    i64* f = out + t.k;
    f[0] = vote, f[1] = spill_n, f[2] = spill, f[3] = miss_n, f[4] = miss;
  }
};

HD inline void sift_min(i32* a, i32 i, i32 n) {
  for (;;) {
    i32 c = 2 * i + 1;
    if (c >= n) return;
    if (c + 1 < n && a[c + 1] < a[c]) ++c;
    if (!(a[c] < a[i])) return;
    const i32 t = a[i];
    a[i] = a[c], a[c] = t, i = c;
  }
}

HD inline void sort_descending(i32* a, i32 n) {  // heapsort on a min-heap: each minimum goes last
  for (i32 i = n / 2 - 1; i >= 0; --i) sift_min(a, i, n);
  for (i32 e = n - 1; e > 0; --e) {
    const i32 t = a[0];
    a[0] = a[e], a[e] = t;
    sift_min(a, 0, e);
  }
}

HD inline double iou(i64 inter, i64 a, i64 b) {
  const i64 uni = a + b - inter;
  return uni ? double(inter) / double(uni) : 1.0;
}

HD inline double sum_value(i128 units) {  // exact below 2^53 units (2^30 px of distance), deterministic above
  const i128 lim = i128(1) << 62;
  if (units < lim && units > -lim) return double(i64(units)) / 8388608.0;
  return (double(i64(units >> 32)) * 4294967296.0 + double(u32(units & 0xFFFFFFFF))) / 8388608.0;
}

// Per object: its signature sorted, then the terms (score_mask's formulas, operation by operation).
struct ObjTerms {
  const Obj* obj;
  const i64* parts;
  const i32* ncomp;
  i32* lists;
  ObjOut* out;
  HD static bool same_sig(const TableView& t, i32 r, const i32* s, i32 m) {
    if (t.sig_start[r + 1] - t.sig_start[r] != m) return false;
    for (i32 j = 0; j < m; ++j)
      if (t.sig[t.sig_start[r] + j] != s[j]) return false;
    return true;
  }
  HD void operator()(i64 o) const {
    const Obj& b = obj[o];
    const TableView& t = b.t;
    const i32 m = ncomp[o], F = t.k + kPart;
    i32* s = lists + b.crop;
    sort_descending(s, m);
    const i64* rows = parts + b.part;
    i64 vote = 0, spill_n = 0, miss_n = 0;
    i128 spill = 0, miss = 0;
    for (i32 y = 0; y < t.H; ++y) {
      const i64* f = rows + i64(y) * F + t.k;
      vote += f[0], spill_n += f[1], spill += f[2], miss_n += f[3], miss += f[4];
    }
    i32 same = 0, voters = 0;
    double best = 0;
    for (i32 r = 0; r < t.k; ++r) {
      if (!t.voter[r]) continue;
      ++voters;
      if (!same_sig(t, r, s, m)) continue;
      i64 inter = 0;
      for (i32 y = 0; y < t.H; ++y) inter += rows[i64(y) * F + r];
      const double v = iou(inter, b.area, t.rater_area[r]);
      if (!same || v > best) best = v;  // Python's max: the first of equal values
      ++same;
    }
    const double d = (sum_value(spill) - 0.5 * double(spill_n)) + (sum_value(miss) - 0.5 * double(miss_n));
    double g = 0.0;
    if (b.area > 0) {
      i64 lo = t.rater_area[0], hi = t.rater_area[0];
      for (i32 r = 1; r < t.k; ++r) lo = mn(lo, t.rater_area[r]), hi = mx(hi, t.rater_area[r]);
      const double up = double(hi) / double(b.area), down = lo > 0 ? double(b.area) / double(lo) : 1.0;
      g = (up < 1.0 ? up : 1.0) * (down < 1.0 ? down : 1.0);
    }
    out[o] = ObjOut{b.area, m, same ? best : iou(vote, b.area, t.n_vote),
                    1.0 / (1.0 + (2.0 * d / double(t.perimeter)) / (t.t5 * t.t5)), g, double(same) / double(voters)};
  }
};

struct SigScan {  // one item
  const i32* ncomp;
  i64 n;
  i64* sig_start;
  HD void operator()(i64) const {
    sig_start[0] = 0;
    for (i64 o = 0; o < n; ++o) sig_start[o + 1] = sig_start[o] + ncomp[o];
  }
};

struct SigCompact {
  const Obj* obj;
  const i32* ncomp;
  const i64* sig_start;
  const i32* lists;
  i32* sig;
  HD void operator()(i64 o) const {
    for (i32 j = 0; j < ncomp[o]; ++j) sig[sig_start[o] + j] = lists[obj[o].crop + j];
  }
};

// Scores the masks at `words` in lane L's checker arena (mask s at word_start[s]) against their tables. The staging
// buffers and the scorer's arena are free (the caller waited for the lane).
void score_words(Lane& L, const std::vector<TableView>& tv, const u32* words, const std::vector<i64>& word_start,
                 Scores& out, Timer& t) {
  const i64 n = i64(tv.size());
  const Stream st = L.stream;
  std::vector<Obj> obj(size_t(n), Obj{});
  std::vector<i64> row_start(size_t(n + 1), 0), crow_start(size_t(n + 1), 0);
  for (i64 o = 0; o < n; ++o) {
    obj[size_t(o)].t = tv[size_t(o)], obj[size_t(o)].words = word_start[size_t(o)];
    row_start[size_t(o + 1)] = row_start[size_t(o)] + tv[size_t(o)].H;
  }
  const i64 rows = row_start[size_t(n)];

  // Boxes: the crop layout.
  Layout a;
  const size_t a_obj = a.put<Obj>(n), a_rows = a.put<i64>(n + 1), a_rbox = a.put<i32>(3 * rows),
               a_box = a.put<i64>(5 * n);
  u8* d = L.score.get(a.size);
  u8* h = L.in.get(a_rbox);
  std::copy_n(obj.data(), n, ptr<Obj>(h, a_obj));
  std::copy_n(row_start.data(), n + 1, ptr<i64>(h, a_rows));
  to_device(d, h, a_rbox, st);
  t.stage("label", st, [&] {  // the crop boxes, which the labeling works on
    for_each(st, rows, RowBox{ptr<Obj>(d, a_obj), ptr<i64>(d, a_rows), n, words, ptr<i32>(d, a_rbox)});
    for_each(st, n, ObjBox{ptr<i64>(d, a_rows), ptr<i32>(d, a_rbox), ptr<i64>(d, a_box)});
  });
  const i64* box = ptr<i64>(L.out.get(sizeof(i64) * size_t(5 * n)), 0);
  to_host(L.out.p, d + a_box, sizeof(i64) * size_t(5 * n), st);
  wait(st);
  i64 parts = 0, crop = 0;
  for (i64 o = 0; o < n; ++o) {
    Obj& b = obj[size_t(o)];
    const i64* x = box + 5 * o;
    b.bx = i32(x[0]), b.by = i32(x[1]), b.bw = i32(x[2]), b.bh = i32(x[3]), b.area = x[4];
    b.part = parts, b.crop = crop;
    parts += i64(b.t.H) * (b.t.k + kPart);
    if (b.area) crop += i64(b.bw + 2) * (b.bh + 2);
    crow_start[size_t(o + 1)] = crow_start[size_t(o)] + (b.area ? b.bh + 2 : 0);
  }
  const i64 crows = crow_start[size_t(n)];
  counters().label_pairs += crows;  // the second label pair relabels every crop row

  // Labeling, row partials, terms.
  Layout l;
  const size_t o_obj = l.put<Obj>(n), o_rows = l.put<i64>(n + 1), o_crows = l.put<i64>(n + 1),
               o_parts = l.put<i64>(parts), o_img = l.put<u8>(crop), o_fg = l.put<i32>(crop), o_bg = l.put<i32>(crop),
               o_aux = l.put<i32>(crop), o_nc = l.put<i32>(n), o_out = l.put<ObjOut>(n), o_ss = l.put<i64>(n + 1);
  d = L.score.get(l.size);
  h = L.in.get(o_parts);
  std::copy_n(obj.data(), n, ptr<Obj>(h, o_obj));
  std::copy_n(row_start.data(), n + 1, ptr<i64>(h, o_rows));
  std::copy_n(crow_start.data(), n + 1, ptr<i64>(h, o_crows));
  to_device(d, h, o_parts, st);
  zero(d + o_nc, sizeof(i32) * size_t(n), st);
  const Obj* dobj = ptr<Obj>(d, o_obj);
  const Crops c{dobj, ptr<i64>(d, o_crows), n, ptr<u8>(d, o_img), ptr<i32>(d, o_fg), ptr<i32>(d, o_bg),
                ptr<i32>(d, o_aux), ptr<i32>(d, o_nc)};
  auto label = [&](int which) {
    t.stage("label", st, [&] {
      for_each(st, crows, LabelInit{c, which});
      for_each(st, crows, LabelMerge{c, which});
      for_each(st, crows, LabelFlatten{c, which});
    });
  };
  auto count = [&](int which) {
    t.stage("count", st, [&] {
      zero(c.aux, sizeof(i32) * size_t(crop), st);
      for_each(st, crows, Count{c, which});
    });
  };
  auto filter = [&](int which) { t.stage("filter", st, [&] { for_each(st, crows, Filter{c, which}); }); };
  t.stage("label", st, [&] { for_each(st, crows, CropFill{c, words}); });
  label(1), count(1), filter(1);
  label(0), count(0), filter(0);
  label(1), label(0);
  t.stage("score", st, [&] {  // holes per component, row partials, the terms and the signatures
    zero(c.aux, sizeof(i32) * size_t(crop), st);
    for_each(st, crows, Holes{c});
    for_each(st, crows, Gather{c});
    for_each(st, rows, ScoreRows{dobj, ptr<i64>(d, o_rows), n, words, ptr<i64>(d, o_parts)});
    for_each(st, n, ObjTerms{dobj, ptr<i64>(d, o_parts), c.ncomp, c.bg, ptr<ObjOut>(d, o_out)});
    for_each(st, 1, SigScan{c.ncomp, n, ptr<i64>(d, o_ss)});
    for_each(st, n, SigCompact{dobj, c.ncomp, ptr<i64>(d, o_ss), c.bg, c.aux});
  });
  const size_t res = sizeof(ObjOut) * size_t(n), res_pad = (res + 255) / 256 * 256;
  u8* ho = L.out.get(res_pad + sizeof(i64) * size_t(n + 1));
  to_host(ho, d + o_out, res, st);
  to_host(ho + res_pad, d + o_ss, sizeof(i64) * size_t(n + 1), st);
  wait(st);
  const ObjOut* r = ptr<ObjOut>(ho, 0);
  for (i64 o = 0; o < n; ++o) {
    out.area.push_back(r[o].area);
    out.iou_tc.push_back(r[o].iou_tc), out.f5.push_back(r[o].f5), out.g_area.push_back(r[o].g_area);
    out.topo_freq.push_back(r[o].topo_freq);
  }
  out.sig_start.assign(ptr<i64>(ho + res_pad, 0), ptr<i64>(ho + res_pad, 0) + n + 1);

  // The signatures.
  const i64 total = out.sig_start.back();
  out.sig.resize(size_t(total));
  to_host(L.out.get(sizeof(i32) * size_t(total)), c.aux, sizeof(i32) * size_t(total), st);
  wait(st);
  std::copy_n(ptr<i32>(L.out.p, 0), total, out.sig.data());
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

Stats stats() {
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  return counters();
}

void reset_stats() {
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Stats& s = counters();
  const Stats kept = s;
  s = Stats{};
  s.arena_bytes = kept.arena_bytes, s.table_bytes = kept.table_bytes;
  s.live_device_bytes = kept.live_device_bytes, s.contexts = kept.contexts;
}

PolygonResult check_polygons(const Polygons& in, Times* times) {
  GC_RANGE("check_polygons");
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
  Timer t{times, {}};
  e.start();
  for (int k = 0; k < used; ++k) {  // every buffer is grown before any work is queued
    Lane& L = e.lanes[k];
    L.dev.get(chunks[k].size), L.in.get(chunks[k].in_bytes), L.out.get(chunks[k].size - chunks[k].defects);
  }
  for (int k = 0; k < used; ++k) {  // lane 1 is packed while lane 0 runs
    pack(in, out.word_start, chunks[k], e.lanes[k].in.p);
    launch(chunks[k], in.G, e.lanes[k], true, t);
  }
  for (int k = 0; k < used; ++k) {
    const PolyChunk& c = chunks[k];
    Lane& L = e.lanes[k];
    wait(L.stream);
    std::copy_n(ptr<i64>(L.out.p, 0), c.n * kNumDefects, out.defects.data() + c.a * kNumDefects);
    std::copy_n(ptr<u32>(L.out.p, c.words - c.defects), c.nwords, out.evenodd.data() + out.word_start[c.a]);
  }
  if (times) t.read();
  return out;
}

MaskResult check_masks(const Masks& in, Times* times) {
  GC_RANGE("check_masks");
  MaskResult out;
  out.ring_start = {0};
  out.pts_start = {0};
  if (in.n == 0) return out;
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  Timer t{times, {}};
  e.start();
  for (i64 m0 = 0; m0 < in.n;) {  // chunks of at most kChunkSlots slots (a larger mask goes alone)
    i64 m1 = m0, slots = 0;
    while (m1 < in.n) {
      const i64 s = 4 * (i64(in.W[m1]) + 1) * (i64(in.H[m1]) + 1);
      if (m1 > m0 && slots + s > kChunkSlots) break;
      slots += s, ++m1;
    }
    check_mask_chunk(in, m0, m1, out, e.lanes[0], t);
    m0 = m1;
  }
  if (times) t.read();
  return out;
}

i64 add_tables(const Tables& t) {
  GC_RANGE("add_tables");
  const i64 wpr = words_per_row(t.W), plane = i64(t.H) * wpr, px = i64(t.H) * t.W, nsig = t.sig_start[t.k];
  Layout l;
  const size_t o_r = l.put<u32>(t.k * plane), o_u = l.put<u32>(plane), o_c = l.put<u32>(plane),
               o_v = l.put<u32>(plane), o_d = l.put<float>(px), o_a = l.put<i64>(t.k), o_ss = l.put<i64>(t.k + 1),
               o_s = l.put<i32>(nsig), o_vt = l.put<u8>(t.k);
  std::vector<u8> h(l.size, 0);
  for (i32 r = 0; r < t.k; ++r) pack_plane(t.raters + r * px, t.H, t.W, ptr<u32>(h.data(), o_r) + r * plane);
  pack_plane(t.uni, t.H, t.W, ptr<u32>(h.data(), o_u));
  pack_plane(t.core, t.H, t.W, ptr<u32>(h.data(), o_c));
  pack_plane(t.vote, t.H, t.W, ptr<u32>(h.data(), o_v));
  std::copy_n(t.dist, px, ptr<float>(h.data(), o_d));
  std::copy_n(t.rater_area, t.k, ptr<i64>(h.data(), o_a));
  std::copy_n(t.sig_start, t.k + 1, ptr<i64>(h.data(), o_ss));
  std::copy_n(t.sig, nsig, ptr<i32>(h.data(), o_s));
  std::copy_n(t.voter, t.k, ptr<u8>(h.data(), o_vt));
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  e.start();
  u8* d = static_cast<u8*>(alloc(true, l.size));
  try {
    to_device(d, h.data(), l.size, e.lanes[0].stream);
    wait(e.lanes[0].stream);
  } catch (...) {
    release(true, d, l.size);
    throw;
  }
  const TableView v{ptr<u32>(d, o_r), ptr<u32>(d, o_u), ptr<u32>(d, o_c), ptr<u32>(d, o_v), ptr<float>(d, o_d),
                    ptr<i64>(d, o_a), ptr<i64>(d, o_ss), ptr<i32>(d, o_s), ptr<u8>(d, o_vt), t.n_vote, t.perimeter,
                    t.t5, t.H, t.W, t.k, t.min_area};
  Registry& r = registry();
  i64 handle = i64(r.block.size());
  if (r.free.empty()) {
    r.block.push_back(d), r.bytes.push_back(l.size), r.view.push_back(v);
  } else {
    handle = r.free.back();
    r.free.pop_back();
    r.block[size_t(handle)] = d, r.bytes[size_t(handle)] = l.size, r.view[size_t(handle)] = v;
  }
  ++r.live;
  counters().table_bytes += i64(l.size);
  return handle;
}

void drop_tables(i64 handle) {
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Registry& r = registry();
  if (handle < 0 || handle >= i64(r.block.size()) || !r.block[size_t(handle)])
    throw std::invalid_argument("unknown tables handle");
  clear_error();
  release(true, r.block[size_t(handle)], r.bytes[size_t(handle)]);
  counters().table_bytes -= i64(r.bytes[size_t(handle)]);
  r.block[size_t(handle)] = nullptr;
  r.free.push_back(handle);
  --r.live;
}

i64 resident() {
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  return registry().live;
}

Scores score_polygons(const Polygons& in, const i64* tables, Times* times) {
  GC_RANGE("score_polygons");
  const i64 n = in.n;
  Scores out;
  out.sig_start = {0};
  std::vector<i64> word_start(size_t(n + 1), 0);
  for (i64 k = 0; k < n; ++k) word_start[size_t(k + 1)] = word_start[size_t(k)] + i64(in.H[k]) * words_per_row(in.W[k]);
  if (n == 0) return out;
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Drain drain{e};
  const std::vector<TableView> tv = views(tables, n, in.W, in.H);
  clear_error();
  Timer t{times, {}};
  e.start();
  Lane& L = e.lanes[0];
  const PolyChunk c = plan(in, word_start, 0, n);
  L.dev.get(c.size), L.in.get(c.in_bytes), L.out.get(c.words - c.defects);
  pack(in, word_start, c, L.in.p);
  launch(c, in.G, L, false, t);
  wait(L.stream);
  out.defects.assign(ptr<i64>(L.out.p, 0), ptr<i64>(L.out.p, 0) + n * kNumDefects);
  score_words(L, tv, ptr<u32>(L.dev.p, c.words), word_start, out, t);
  if (times) t.read();
  return out;
}

Scores score_masks(const Masks& in, const i64* tables, Times* times) {
  GC_RANGE("score_masks");
  const i64 n = in.n;
  Scores out;
  out.sig_start = {0};
  std::vector<i64> word_start(size_t(n + 1), 0);
  for (i64 k = 0; k < n; ++k) word_start[size_t(k + 1)] = word_start[size_t(k)] + i64(in.H[k]) * words_per_row(in.W[k]);
  if (n == 0) return out;
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Drain drain{e};
  const std::vector<TableView> tv = views(tables, n, in.W, in.H);
  clear_error();
  Timer t{times, {}};
  e.start();
  Lane& L = e.lanes[0];
  const size_t bytes = sizeof(u32) * size_t(word_start[size_t(n)]);
  u8* h = L.in.get(bytes);
  std::fill_n(h, bytes, u8(0));
  for (i64 k = 0; k < n; ++k)
    pack_plane(in.bits + in.start[k], in.H[k], in.W[k], ptr<u32>(h, 0) + word_start[size_t(k)]);
  u8* d = L.dev.get(bytes);
  to_device(d, h, bytes, L.stream);
  wait(L.stream);
  out.defects.assign(size_t(n * kNumDefects), 0);
  score_words(L, tv, ptr<u32>(d, 0), word_start, out, t);
  if (times) t.read();
  return out;
}

}  // namespace gpucheck
