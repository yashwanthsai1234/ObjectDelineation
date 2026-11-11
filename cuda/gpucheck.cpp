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
// Masks (objdel.checker.polygonize's rings and bitquad's counts, the K1-draft design; "masks" section below):
//  * bit-packed rows, runs per row, union-find on runs (foreground 8-, background 4-connected), roots = first run;
//  * one ring per foreground component (outer) and per bounded background component (hole): its order key (first east
//    edge in raster order) and canonical start in closed form from the labels, its rank by a popcount scan (no sort);
//  * crack tracing from seeds over boundary edges only, corners written in place: work ~ words + runs + edges.
//
// Host side (issues CK-1 to CK-3 of experiments/phase2/perf/issues.md): per lane one stream, device arenas and two
// staging buffers (pinned host memory). Device memory (arenas, tables) comes from the engine's own stream-ordered pool
// (D5c: allocated and freed in stream order, no device-wide sync; the process's default pool is never touched); an
// arena grows to the next power of two and is never freed per call. The staging buffers are sized once, in the first
// call of the process, to the bounds of a 256-answer score_polygons call (stage bounds below), so the RL step never
// grows them (a growth synchronizes the device; it is counted as pinned_grows and logged). A polygon call is: pack
// the inputs into the pinned buffer, one copy in, the kernels on the lane's stream, one copy out, one stream
// synchronize. A large polygon batch is cut into two halves on two lanes, so the second half's copy and the host's
// packing overlap the first half's kernels (C7). One mutex serializes the calls.
//
// References: the headers of archive/cpu_checker/cpp/polycore/polycore.hpp and src/objdel/checker/polygonize.py; CUDA C++ Programming
// Guide (https://docs.nvidia.com/cuda/cuda-c-programming-guide/): __int128 in device code (CUDA >= 11.5, "Supported
// Types"), warp vote and shuffle functions (__ballot_sync, __shfl_sync: every lane named in the mask must execute
// them), __syncwarp (memory ordering among the warp's lanes), page-locked host memory and asynchronous concurrent
// execution (cudaMemcpyAsync overlaps kernels only from pinned memory, on a non-default stream).
#include "gpucheck.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>

#ifdef __CUDACC__
#include <cuda_runtime.h>

#include <cub/device/device_scan.cuh>
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

// Test-only compile definitions (CMake OBJDEL_TEST_HOOKS, set by tests/variants.py only; test_hooks() names them):
// K1_ORDER 1/2 runs the CPU build's items reversed / permuted, GC_POISON_SCAN fills a scan's scratch with 0xFF,
// GC_ALIAS_SCAN_FOR_TEST hands one scan a live block, K1_MUTATE_FG4 makes the polygonizer foreground 4-connected.
#ifndef K1_ORDER
#define K1_ORDER 0
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

// The engine's device memory pool (made by Engine::start): its own pool on the current device, never set as the
// device's pool (cudaDeviceSetMemPool would change torch's and vLLM's allocations too), release threshold UINT64_MAX
// so a sync never shrinks it (CUDA programming guide, stream-ordered memory allocator; cudaMemPoolCreate,
// cudaMallocFromPoolAsync, cudaMemPoolTrimTo in the runtime API reference).
cudaMemPool_t& pool() {
  static cudaMemPool_t p = nullptr;
  return p;
}

cudaMemPool_t new_pool() {
  int dev = 0;
  cuda_check(cudaGetDevice(&dev));
  cudaMemPoolProps props{};
  props.allocType = cudaMemAllocationTypePinned;
  props.location.type = cudaMemLocationTypeDevice;
  props.location.id = dev;
  cudaMemPool_t p = nullptr;
  cuda_check(cudaMemPoolCreate(&p, &props));
  std::uint64_t keep = UINT64_MAX;
  cuda_check(cudaMemPoolSetAttribute(p, cudaMemPoolAttrReleaseThreshold, &keep));
  return p;
}

i64 reserved(cudaMemPool_t p) {
  std::uint64_t v = 0;
  cuda_check(cudaMemPoolGetAttribute(p, cudaMemPoolAttrReservedMemCurrent, &v));
  return i64(v);
}

// Device memory: from the engine's pool in stream order on `st`. Host staging: cudaMallocHost / cudaFreeHost, which
// synchronize the device (counted as sync_alloc / sync_free).
void* alloc(bool device, size_t bytes, Stream st) {
  void* p = nullptr;
  if (device) {
    cuda_check(cudaMallocFromPoolAsync(&p, bytes, pool(), st));
    ++counters().pool_allocs, counters().live_device_bytes += i64(bytes);
  } else {
    cuda_check(cudaMallocHost(&p, bytes));
    ++counters().sync_alloc;
  }
  return p;
}

void release(bool device, void* p, size_t bytes, Stream st) {
  if (device) {
    ++counters().pool_frees, counters().live_device_bytes -= i64(bytes);
    cuda_check(cudaFreeAsync(p, st));
  } else {
    ++counters().sync_free;
    cuda_check(cudaFreeHost(p));
  }
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

void copy_device(void* dst, const void* src, size_t bytes, Stream st) {
  if (bytes) cuda_check(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToDevice, st));
}

// Device -> pageable host memory after a wait (the device is idle), queued on the lane stream `st`: the pool buffer is
// allocated and freed in stream order on `st` (cudaMallocFromPoolAsync / cudaFreeAsync on a non-blocking stream), so
// a default-stream cudaMemcpy has no order with them (box_a memcheck). cudaMemcpyAsync to pageable memory only "might"
// return after the copy (runtime API, API synchronization behavior); the stream sync makes it complete before `dst`
// is read. Not a counted wait.
void read_back(void* dst, const void* src, size_t bytes, Stream st) {
  if (!bytes) return;
  cuda_check(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, st));
  cuda_check(cudaStreamSynchronize(st));
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
// Counted as the CUDA build counts them: device memory as pool allocations, host staging as synchronizing ones.
void* alloc(bool device, size_t bytes, Stream) {
  void* p = ::operator new(bytes);
  ++(device ? counters().pool_allocs : counters().sync_alloc);
  if (device) counters().live_device_bytes += i64(bytes);
  return p;
}
void release(bool device, void* p, size_t bytes, Stream) {
  ++(device ? counters().pool_frees : counters().sync_free);
  if (device) counters().live_device_bytes -= i64(bytes);
  ::operator delete(p);
}
void to_device(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void to_host(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void zero(void* p, size_t bytes, Stream) { if (bytes) std::memset(p, 0, bytes); }
void copy_device(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void read_back(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void wait(Stream) { ++counters().waits; }
void drain_sync(Stream) { ++counters().drain_syncs; }
void clear_error() {}

// Stage marks (Timer): the host clock, since every step has run when the next one is queued.
using Mark = std::chrono::steady_clock::time_point;
Mark mark(Stream) { return std::chrono::steady_clock::now(); }
double elapsed_ms(Mark a, Mark b) { return std::chrono::duration<double, std::milli>(b - a).count(); }
void unmark(Mark) {}

// K1_ORDER stands in for the GPU's arbitrary schedule: 0 forward, 1 reverse, 2 the permutation k -> k p mod n
// (p = 2^31 - 1 is prime, so coprime with any n < p). Every step must give the same result in every order.
template <class F>
void for_each(Stream, i64 n, const F& f) {
  for (i64 k = 0; k < n; ++k) f(K1_ORDER == 0 ? k : K1_ORDER == 1 ? n - 1 - k : (k * 2147483647LL) % n);
}
#endif

// The entry and size of the running call (set by Engine::start), named by the staging growth log.
struct Call {
  const char* entry = "";
  i64 n = 0;
};
Call& call() {
  static Call c;
  return c;
}

// Grow-only buffer: device memory (stream-ordered on `st`, the lane's stream), or host staging memory (pinned on the
// GPU build). get(bytes) keeps the buffer when it is large enough, else replaces it by one of the next power of two up
// to 64 MiB, exact above (the polygon path stays in a few MiB; one large check_masks call must not keep twice its
// need). Never freed per call. A staging buffer is sized once by Engine::start (resize), so its growth in get is the
// event pinned_grows, logged with the call that needed it.
template <bool kDevice>
struct Block {
  u8* p = nullptr;
  size_t cap = 0;
  Stream st{};
  void resize(size_t c) {
    void* old = p;
    const size_t old_cap = cap;
    p = nullptr, cap = 0;
    if (kDevice) counters().arena_bytes -= i64(old_cap);
    if (old) release(kDevice, old, old_cap, st);
    p = static_cast<u8*>(alloc(kDevice, c, st));
    cap = c;
    if (kDevice) counters().arena_bytes += i64(cap);
  }
  u8* get(size_t bytes) {
    if (bytes <= cap) return p;
    size_t c = 4096;
    while (c < bytes && c < (size_t(64) << 20)) c <<= 1;
    if (c < bytes) c = bytes;
    if (!kDevice) {
      ++counters().pinned_grows;
      std::fprintf(stderr, "gpucheck: pinned staging grows from %zu to %zu bytes in %s (n = %lld, %zu bytes asked)\n",
                   cap, c, call().entry, static_cast<long long>(call().n), bytes);
    }
    resize(c);
    return p;
  }
};

struct Lane {  // one stream with its device arenas (dev: checker, score: scorer; the polygonizer's layouts A and B),
  Stream stream{};  // CUB's scan scratch (scan, never live data) and its staging buffers (in: host -> device, out:
  Block<true> dev, score, scan;  // device -> host)
  Block<false> in, out;
};

size_t stage_in_bound();  // the staging bounds (after the polygon layout, below)
size_t stage_out_bound();

struct Engine {
  std::mutex mu;
  bool ready = false;
  Lane lanes[2];
  // Called with mu held by every entry before it stages a byte: names the call; the first call of the process makes
  // the streams and the pool and sizes every lane's staging buffers to their bounds.
  void start(const char* entry, i64 n) {
    call() = Call{entry, n};
    if (ready) return;
    for (Lane& l : lanes) l.stream = l.dev.st = l.score.st = l.scan.st = new_stream();
#ifdef __CUDACC__
    if (!pool()) pool() = new_pool();
#endif
    for (Lane& l : lanes) l.in.resize(stage_in_bound()), l.out.resize(stage_out_bound());
    ready = true;
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

HD inline void atomic_or(u32* p, u32 v) {
#ifdef __CUDA_ARCH__
  atomicOr(p, v);
#else
  *p |= v;
#endif
}

HD inline void set_flag(i32* p) {  // an error flag any thread may set
#ifdef __CUDA_ARCH__
  atomicExch(p, 1);
#else
  *p = 1;
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
template <class T>
HD inline i64 segment(const T* start, i64 n, i64 i) {
  i64 lo = 0, hi = n;
  while (hi - lo > 1) {
    const i64 mid = (lo + hi) / 2;
    if (start[mid] <= i) lo = mid;
    else hi = mid;
  }
  return lo;
}

HD inline i64 iabs(i64 v) { return v < 0 ? -v : v; }

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

// Word offsets of n bit-packed frames: frame k is words [ws[k], ws[k + 1]), H[k] rows of words_per_row(W[k]).
std::vector<i64> word_starts(i64 n, const i32* W, const i32* H) {
  std::vector<i64> ws(size_t(n + 1), 0);
  for (i64 k = 0; k < n; ++k) ws[size_t(k + 1)] = ws[size_t(k)] + i64(H[k]) * words_per_row(W[k]);
  return ws;
}

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

// ---------------------------------------------------------------- warp helpers (K2-draft k2score.inc)
// for_each_warp runs one item per warp of 32 lanes on the GPU and one item per "warp of one lane" on the CPU build, so
// lane loops (for i = lane; i < n; i += size) become plain loops there and the collectives the identity. Every lane of
// the warp must reach each collective (warp-uniform control flow up to the call; CUDA C++ guide, warp intrinsics).
struct Warp {
  int lane, size;  // GPU: lane in [0, 32), size 32. CPU build: lane 0, size 1
};

#ifdef __CUDACC__
template <class F>
__global__ void warp_kernel(i64 n, F f) {
  const i64 i = (i64(blockIdx.x) * blockDim.x + threadIdx.x) >> 5;  // warp-uniform: whole warps return together
  if (i < n) f(i, Warp{int(threadIdx.x & 31), 32});
}

template <class F>
void for_each_warp(Stream st, i64 n, const F& f) {
  constexpr int kThreads = 128;  // 4 warps per block
  if (n <= 0) return;
  warp_kernel<<<unsigned((n * 32 + kThreads - 1) / kThreads), kThreads, 0, st>>>(n, f);
  cuda_check(cudaGetLastError());
}
#else
template <class F>
void for_each_warp(Stream, i64 n, const F& f) {
  for (i64 k = 0; k < n; ++k) f(K1_ORDER == 0 ? k : K1_ORDER == 1 ? n - 1 - k : (k * 2147483647LL) % n, Warp{0, 1});
}
#endif

HD inline u32 wballot(bool b) {
#ifdef __CUDA_ARCH__
  return __ballot_sync(0xFFFFFFFFu, b);
#else
  return b ? 1u : 0u;
#endif
}
HD inline bool wany(bool b) { return wballot(b) != 0; }
HD inline bool wodd(bool b) { return (popc(wballot(b)) & 1) != 0; }  // parity of the lanes' flags
HD inline i64 wsum(i64 v) {  // integers: the same value in any order (__reduce_add_sync takes 32-bit values only)
#ifdef __CUDA_ARCH__
  for (int o = 16; o > 0; o >>= 1) v += __shfl_xor_sync(0xFFFFFFFFu, v, o);
#endif
  return v;
}
HD inline void wsync() {  // orders the lanes' global-memory writes before the next phase's reads
#ifdef __CUDA_ARCH__
  __syncwarp();
#endif
}

// Winding and containment of one shape, one warp per shape (polycore detail::validate_geometry, second half; K2-draft
// ShapeTotalsWarp): the 13 totals, written whole. Phase 1, one lane per ring: wrong winding. Phase 2, the whole warp
// per containment test: many witnesses take one lane each (the serial locate), 32 at a time, leaving at the first
// chunk with a hit; a single witness splits the ring's edges over the lanes. Phase 3: totals by warp sums. Every
// result is an integer or a flag, so the lane split and the order do not change it.
struct ShapeTotalsWarp {
  Shapes s;
  Rings r;
  const u8* validate_wide;
  i64* defects;

  // Point in doubled coordinates against clean ring g: -1 strictly outside, 0 on the boundary, 1 inside. The ring's
  // edges are split over w's lanes; kWarp combines their flags (on = any lane, inside = parity of the crossings), else
  // one lane (w = {0, 1}) runs every edge.
  template <class T, bool kWarp>
  HD int locate(i64 g, i64 qx, i64 qy, Warp w) const {
    const i32* box = r.box + 4 * g;
    if (qx < 2 * i64(box[0]) || qx > 2 * i64(box[1]) || qy < 2 * i64(box[2]) || qy > 2 * i64(box[3])) return -1;
    const i32* c = r.cxy + 2 * s.ring_start[g];
    const i32 n = r.len[g];
    bool on = false, in = false;
    for (i32 k = w.lane; k < n; k += w.size) {
      const i32 k1 = k + 1 == n ? 0 : k + 1;
      const i64 ax = 2 * i64(c[2 * k]), ay = 2 * i64(c[2 * k + 1]), bx = 2 * i64(c[2 * k1]), by = 2 * i64(c[2 * k1 + 1]);
      if (in_box(ax, ay, bx, by, qx, qy) && orient<T>(ax, ay, bx, by, qx, qy) == 0) {
        on = true;
      } else if ((ay <= qy && qy < by) || (by <= qy && qy < ay)) {
        const T num = T(ax - qx) * T(by - ay) + T(qy - ay) * T(bx - ax);
        if (by > ay ? num > 0 : num < 0) in = !in;
      }
    }
    if (kWarp) on = wany(on), in = wodd(in);
    return on ? 0 : in ? 1 : -1;
  }

  // Region of outer ring o (local) with its holes, even-odd over all of them.
  template <class T, bool kWarp>
  HD int locate_region(i64 r0, i32 R, i32 o, i64 qx, i64 qy, Warp w) const {
    const int loc = locate<T, kWarp>(r0 + o, qx, qy, w);
    if (loc == 0) return 0;
    bool in = loc > 0;
    for (i32 h = 0; h < R; ++h) {
      if (s.roles[r0 + h] != o || r.fewer[r0 + h]) continue;
      const int lh = locate<T, kWarp>(r0 + h, qx, qy, w);
      if (lh == 0) return 0;
      if (lh > 0) in = !in;
    }
    return in ? 1 : -1;
  }

  // Does a witness of ring a (vertex or edge midpoint) locate to `want` against ring b (or b's region)?
  // single: a neither crosses nor touches that boundary, so its first vertex decides.
  template <class T>
  HD bool any_witness(i64 r0, i32 R, i32 a, int want, bool single, i32 b, bool region, Warp w) const {
    const i32* c = r.cxy + 2 * s.ring_start[r0 + a];
    const i32 n = r.len[r0 + a];
    if (single) {
      if (n <= 0) return false;
      const i64 x = 2 * i64(c[0]), y = 2 * i64(c[1]);
      return (region ? locate_region<T, true>(r0, R, b, x, y, w) : locate<T, true>(r0 + b, x, y, w)) == want;
    }
    const Warp one{0, 1};
    for (i32 base = 0; base < 2 * n; base += w.size) {  // warp-uniform trip count
      const i32 i = base + w.lane;
      bool hit = false;
      if (i < 2 * n) {
        const i32 k = i >> 1, k1 = k + 1 == n ? 0 : k + 1;
        const i64 x = (i & 1) ? i64(c[2 * k]) + c[2 * k1] : 2 * i64(c[2 * k]);
        const i64 y = (i & 1) ? i64(c[2 * k + 1]) + c[2 * k1 + 1] : 2 * i64(c[2 * k + 1]);
        hit = (region ? locate_region<T, false>(r0, R, b, x, y, one) : locate<T, false>(r0 + b, x, y, one)) == want;
      }
      if (wany(hit)) return true;
    }
    return false;
  }

  template <class T>
  HD void run(i64 sh, Warp w) const {
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
    for (i32 a = w.lane; a < R; a += w.size) {  // phase 1
      if (fewer(a)) continue;
      const i64 a2 = r.area2[r0 + a];
      if ((roles[a] < 0 && a2 < 0) || (roles[a] >= 0 && a2 > 0)) at(a, kWrongWinding) = 1;
    }
    i64 nested_h = 0, nested_o = 0;  // phase 2: uniform loops, lane 0 writes
    for (i32 h = 0; h < R; ++h) {
      const i32 p = roles[h];
      if (p < 0 || fewer(h) || fewer(p)) continue;
      if (any_witness<T>(r0, R, h, -1, apart(h, p), p, false, w) && w.lane == 0) at(h, kHoleOutsideOuter) = 1;
    }
    for (i32 a = 0; a < R; ++a) {
      if (fewer(a)) continue;
      for (i32 b = a + 1; b < R; ++b) {
        if (fewer(b)) continue;
        if (roles[a] >= 0 && roles[a] == roles[b]) {  // sibling holes
          const bool single = apart(a, b);
          if (any_witness<T>(r0, R, b, 1, single, a, false, w) || any_witness<T>(r0, R, a, 1, single, b, false, w))
            ++nested_h;
        } else if (roles[a] < 0 && roles[b] < 0) {  // two outers
          if (any_witness<T>(r0, R, b, 1, apart_region(b, a), a, true, w) ||
              any_witness<T>(r0, R, a, 1, apart_region(a, b), b, true, w))
            ++nested_o;
        }
      }
    }
    wsync();
    i64* tot = defects + sh * kNumDefects;  // phase 3
    for (int d = 0; d < kNumDefects; ++d) {
      i64 v = 0;
      if (d == kRingRingCrossing || d == kRingRingTouching) {
        const i32* m = d == kRingRingCrossing ? pcross : ptouch;
        for (i32 q = w.lane; q < R * R; q += w.size)
          if (q / R < q % R) v += m[(q / R) * kMaxRings + q % R];
      } else if (d == kNestedHoles || d == kNestedOuters) {
        v = w.lane == 0 ? (d == kNestedHoles ? nested_h : nested_o) : 0;
      } else {
        for (i32 a = w.lane; a < R; a += w.size) v += at(a, d);
      }
      v = wsum(v);
      if (w.lane == 0) tot[d] = v;
    }
  }

  HD void operator()(i64 sh, Warp w) const {
    if (validate_wide[sh]) run<i128>(sh, w);
    else run<i64>(sh, w);
  }
};

// ===================================================================== masks
// The K1-draft polygonizer (docs/implementation2.0/kernels/K1, DESIGN.md: proofs, complexity, decision record) on a
// lane. Conventions of objdel.checker.polygonize and bitquad: vertex (x, y), 0 <= x <= W, 0 <= y <= H; edge = (vertex,
// direction d), d 0 east, 1 south, 2 west, 3 north, key (y (W + 1) + x) 4 + d; foreground lies right of every edge; the
// bit-quad at vertex (x, y) is TL = pixel (x-1, y-1) = 1, TR = (x, y-1) = 2, BL = (x-1, y) = 4, BR = (x, y) = 8; at a
// saddle the walk turns left (8-connected foreground, 4-connected background). Steps of one chunk of masks:
//   Pack      bytes -> bit words (pixel x = bit x % 32 of word x / 32), one item per word (a warp ballot on the GPU).
//   Quads     bit-quad counts n1, n2, nD, n3, n4 by bit slicing on words, one item per vertex row.
//   RunCount  foreground / background runs per row; scan -> run offsets. [host wait 1: the run totals size the rest]
//   RunEmit   runs (x0, x1) per row, in raster order; background run 0 of every mask is its OUTSIDE sentinel.
//   Merge     union-find on runs: foreground rows joined if 8-adjacent, background if 4-adjacent, border background
//             to the sentinel; roots are the minimum index = the first run in raster order. Flatten: run -> root.
//   HoleKey   per bounded background component B (= one hole ring) the raster index of its ring's first east edge.
//   RingMark  order bitmap (one bit per ring at its order pixel), start bitmap (at its start pixel), Betti numbers.
//   scans     bitmap popcounts -> ring index (polygonize's order, no sort), start rank, outer-start rank (parts).
//   RingInit  each ring's start edge in closed form, and its part (raster rank of its foreground component).
//   Seeds     the ring starts plus, on every seed_rows-th row, both vertical sides of each foreground run.
//   Trace     pass 1: each seed walks boundary edges to the next seed (corners, area, next seed); RingWalk: each
//             segment's offset in its ring; scan -> output offsets; pass 2: each seed writes its corners in place.
//   [host wait 2: the ring and corner totals; host wait 3: the outputs]
// Every value is an integer, and cross-thread writes are commuting integer atomics, union-find links (root = minimum)
// or an item's own slots, so the result does not depend on the schedule (the K1_ORDER test builds run the CPU build's
// items reversed and permuted).

// Exclusive prefix sum of n items (out[k] = in[0] + ... + in[k-1]); callers pass one trailing zero to get the total.
// `tmp` is CUB's scratch: it must never be a block that holds live data (the lane's `scan` arena only).
void scan(Stream st, const i32* in, i32* out, i64 n, Block<true>& tmp) {
  if (n <= 0) return;
  if (n > INT_MAX) throw std::logic_error("scan: n > INT_MAX");  // the i32 sums would wrap
#ifdef __CUDACC__
  size_t bytes = 0;
  cuda_check(cub::DeviceScan::ExclusiveSum(nullptr, bytes, in, out, n, st));
  void* t = tmp.get(bytes > 0 ? bytes : 1);
  cuda_check(cub::DeviceScan::ExclusiveSum(t, bytes, in, out, n, st));
#else
  (void)st;
#ifdef GC_POISON_SCAN  // test build: scratch is overwritten, as CUB would, so a live block passed as `tmp` breaks
  std::memset(tmp.get(4 * size_t(n + 1)), 0xFF, tmp.cap);
#else
  (void)tmp;
#endif
  i32 s = 0;
  for (i64 k = 0; k < n; ++k) {
    const i32 v = in[k];
    out[k] = s;
    s += v;
  }
#endif
}

HD inline i32 step_x(int d) { return d == 0 ? 1 : (d == 2 ? -1 : 0); }
HD inline i32 step_y(int d) { return d == 1 ? 1 : (d == 3 ? -1 : 0); }
HD inline bool has_edge(int code, int d) {
  switch (d) {
    case 0: return (code & 8) && !(code & 2);   // top side of BR
    case 1: return (code & 4) && !(code & 8);   // right side of BL
    case 2: return (code & 1) && !(code & 4);   // bottom side of TL
    default: return (code & 2) && !(code & 1);  // left side of TR
  }
}
HD inline i64 area_term(i32 x, i32 y, int d) {  // x y1 - x1 y of the unit edge from (x, y) in direction d
  return d == 0 ? -i64(y) : d == 1 ? i64(x) : d == 2 ? i64(y) : -i64(x);
}

struct Geo {  // one chunk of masks
  const i32 *W, *H, *wpr;        // per mask: width, height, words per row
  const i32 *ws, *rs, *vs, *ss;  // prefix (n + 1): words, rows, vertex rows (H + 1), seed-map entries
  const i32* rowmask;            // per row: its mask
  i32 n, K;                      // masks, seed rows (0 = off)

  HD u32 word(const u32* bits, i32 m, i32 y, i32 j) const {  // word j of row y, 0 outside the mask
    return (y >= 0 && y < H[m] && j >= 0 && j < wpr[m]) ? bits[ws[m] + y * wpr[m] + j] : 0u;
  }
  HD bool px(const u32* bits, i32 m, i32 x, i32 y) const {
    if (x < 0 || y < 0 || x >= W[m] || y >= H[m]) return false;
    return (bits[ws[m] + y * wpr[m] + (x >> 5)] >> (x & 31)) & 1u;
  }
  HD int quad(const u32* bits, i32 m, i32 x, i32 y) const {
    return px(bits, m, x - 1, y - 1) | px(bits, m, x, y - 1) << 1 | px(bits, m, x - 1, y) << 2 | px(bits, m, x, y) << 3;
  }
  // Direction of the edge after the one arriving at vertex (x, y) in direction d: the left turn if it exists, else
  // straight, else right (polygonize's link).
  HD int succ_dir(const u32* bits, i32 m, i32 x, i32 y, int d) const {
    const int c = quad(bits, m, x, y);
#ifdef K1_MUTATE_FG4  // test build (D-G16): the right turn first, foreground 4-connected
    const int turns[3] = {(d + 1) & 3, d, (d + 3) & 3};
#else
    const int turns[3] = {(d + 3) & 3, d, (d + 1) & 3};
#endif
    for (int t : turns)
      if (has_edge(c, t)) return t;
    return -1;
  }
  // Direction of the edge arriving at vertex (x, y) whose successor leaves in direction t (the successor map is a
  // permutation of the edges, so exactly one exists).
  HD int pred_dir(const u32* bits, i32 m, i32 x, i32 y, int t) const {
    for (int d = 0; d < 4; ++d)
      if (has_edge(quad(bits, m, x - step_x(d), y - step_y(d)), d) && succ_dir(bits, m, x, y, d) == t) return d;
    return -1;
  }
  HD u32 valid(i32 m, i32 j) const {  // bits of word j inside the row
    const i32 r = W[m] - 32 * j;
    return r >= 32 ? ~0u : (r <= 0 ? 0u : ((1u << r) - 1u));
  }
  HD i32 bitword(i32 m, i32 x, i32 y) const { return ws[m] + y * wpr[m] + (x >> 5); }
};

struct Runs {
  const i32 *rf, *rb;                           // per row (rows + 1): first foreground / background run
  i32 *fx0, *fx1, *frow, *fpar;                 // foreground runs
  i32 *bx0, *bx1, *brow, *bpar, *bleft, *hkey;  // background runs; bx0 = -1 marks a mask's OUTSIDE sentinel
};

// ---------------------------------------------------------------- words, quads, runs
struct Pack {
  Geo g;
  const u8* bytes;
  const i64* start;  // chunk-local byte offset of each mask
  u32* words;
  HD void operator()(i64 w) const {
    const i32 m = i32(segment(g.ws, g.n, w));
    const i32 lw = i32(w - g.ws[m]), wpr = g.wpr[m], W = g.W[m];
    const i32 y = lw / wpr, x0 = 32 * (lw - y * wpr);
    const u8* row = bytes + start[m] + i64(y) * W;
    u32 v = 0;
    for (i32 b = 0; b < mn(W - x0, 32); ++b) v |= u32(row[x0 + b] != 0) << b;
    words[w] = v;
  }
};

#ifdef __CUDACC__
// One warp per word: lane l reads pixel 32 j + l (coalesced bytes) and __ballot_sync builds the word (as the
// raster). The early exit is uniform per warp (blockDim is a multiple of 32 and all 32 lanes of a warp share w), so
// every lane named in the ballot mask executes it.
__global__ void pack_kernel(Pack p, i64 nw) {
  const i64 t = i64(blockIdx.x) * blockDim.x + threadIdx.x, w = t >> 5;
  const int lane = int(threadIdx.x & 31u);
  if (w >= nw) return;
  const Geo& g = p.g;
  const i32 m = i32(segment(g.ws, g.n, w));
  const i32 lw = i32(w - g.ws[m]), wpr = g.wpr[m], W = g.W[m];
  const i32 y = lw / wpr, x = 32 * (lw - y * wpr) + lane;
  const bool v = x < W && p.bytes[p.start[m] + i64(y) * W + x] != 0;
  const u32 word = __ballot_sync(0xffffffffu, v);
  if (lane == 0) p.words[w] = word;
}
void pack_words(Stream st, const Pack& p, i64 nw) {
  constexpr int kThreads = 256;
  if (nw <= 0) return;
  pack_kernel<<<unsigned((nw * 32 + kThreads - 1) / kThreads), kThreads, 0, st>>>(p, nw);
  cuda_check(cudaGetLastError());
}
#else
void pack_words(Stream st, const Pack& p, i64 nw) { for_each(st, nw, p); }
#endif

struct Quads {  // bit slicing over 32 vertex columns at once
  Geo g;
  const u32* words;
  i64* quads;
  HD void operator()(i64 i) const {
    const i32 m = i32(segment(g.vs, g.n, i));
    const i32 y = i32(i - g.vs[m]), W = g.W[m];
    i64 n1 = 0, n2 = 0, nd = 0, n3 = 0, n4 = 0;
    for (i32 j = 0; j < (W + 1 + 31) / 32; ++j) {
      const u32 B = g.word(words, m, y - 1, j), D = g.word(words, m, y, j);  // TR, BR: pixel x
      const u32 A = (B << 1) | (g.word(words, m, y - 1, j - 1) >> 31);       // TL: pixel x - 1
      const u32 C = (D << 1) | (g.word(words, m, y, j - 1) >> 31);           // BL
      const i32 r = W + 1 - 32 * j;                                         // vertex columns x <= W
      const u32 vm = r >= 32 ? ~0u : ((1u << r) - 1u);
      const u32 e1 = (A & ~B & ~C & ~D) | (~A & B & ~C & ~D) | (~A & ~B & C & ~D) | (~A & ~B & ~C & D);
      const u32 e3 = (~A & B & C & D) | (A & ~B & C & D) | (A & B & ~C & D) | (A & B & C & ~D);
      const u32 e4 = A & B & C & D;
      const u32 dg = (A & D & ~B & ~C) | (B & C & ~A & ~D);  // codes 9 and 6
      const u32 e2 = (A | B | C | D) & ~e1 & ~e3 & ~e4;
      n1 += popc(e1 & vm), n2 += popc(e2 & vm & ~dg), nd += popc(dg & vm), n3 += popc(e3 & vm), n4 += popc(e4 & vm);
    }
    i64* q = quads + 5 * i64(m);
    atomic_add(q + 0, n1), atomic_add(q + 1, n2), atomic_add(q + 2, nd), atomic_add(q + 3, n3), atomic_add(q + 4, n4);
  }
};

struct RunCount {
  Geo g;
  const u32* words;
  i32 *cf, *cb;
  HD void operator()(i64 r) const {
    const i32 m = g.rowmask[r], y = i32(r - g.rs[m]);
    u32 pf = 0, pb = 0;
    i32 nf = 0, nbg = 0;
    for (i32 j = 0; j < g.wpr[m]; ++j) {
      const u32 w = g.word(words, m, y, j), nw = ~w & g.valid(m, j);
      nf += popc(w & ~((w << 1) | pf));
      nbg += popc(nw & ~((nw << 1) | pb));
      pf = w >> 31, pb = nw >> 31;
    }
    cf[r] = nf;
    cb[r] = nbg + (y == 0 ? 1 : 0);  // the mask's OUTSIDE sentinel takes the first background slot of row 0
  }
};

struct RunEmit {
  Geo g;
  const u32* words;
  Runs R;
  HD void operator()(i64 r) const {
    const i32 m = g.rowmask[r], y = i32(r - g.rs[m]), wpr = g.wpr[m];
    i32 bi = R.rb[r];
    if (y == 0) {
      R.bx0[bi] = -1, R.bx1[bi] = -1, R.brow[bi] = i32(r), R.bpar[bi] = bi, R.bleft[bi] = -1, R.hkey[bi] = INT_MAX;
      ++bi;
    }
    const i32 f0 = R.rf[r], b0 = bi;
    i32 fs = 0, fe = 0, bs = 0, be = 0;
    u32 pf = 0, pb = 0;
    for (i32 j = 0; j < wpr; ++j) {
      const u32 w = g.word(words, m, y, j), nw = ~w & g.valid(m, j);
      const u32 wn = g.word(words, m, y, j + 1), nwn = j + 1 < wpr ? (~wn & g.valid(m, j + 1)) : 0u;
      u32 s1 = w & ~((w << 1) | pf), e1 = w & ~((w >> 1) | ((wn & 1u) << 31));
      u32 s0 = nw & ~((nw << 1) | pb), e0 = nw & ~((nw >> 1) | ((nwn & 1u) << 31));
      for (; s1; s1 &= s1 - 1) {
        const i32 k = f0 + fs++;
        R.fx0[k] = 32 * j + ctz(s1), R.frow[k] = i32(r), R.fpar[k] = k;
      }
      for (; e1; e1 &= e1 - 1) R.fx1[f0 + fe++] = 32 * j + ctz(e1);
      for (; s0; s0 &= s0 - 1) {
        const i32 k = b0 + bs++;
        R.bx0[k] = 32 * j + ctz(s0), R.brow[k] = i32(r), R.bpar[k] = k, R.hkey[k] = INT_MAX;
      }
      for (; e0; e0 &= e0 - 1) R.bx1[b0 + be++] = 32 * j + ctz(e0);
      pf = w >> 31, pb = nw >> 31;
    }
    // The foreground run left of each background run (the one ending at x0 - 1): runs alternate along the row.
    i32 a = 0;
    for (i32 k = 0; k < bs; ++k) {
      const i32 x0 = R.bx0[b0 + k];
      while (a < fs && R.fx1[f0 + a] < x0 - 1) ++a;
      R.bleft[b0 + k] = (x0 > 0 && a < fs && R.fx1[f0 + a] == x0 - 1) ? f0 + a : -1;
    }
  }
};

// Union of runs of row y with row y - 1, and of border background runs with the sentinel. Two-pointer merge of two
// sorted lists of disjoint intervals: when two runs end at different x, the one ending first meets no later run of
// the other row (its next run starts >= 2 past the end); equal ends advance both.
struct Merge {
  Geo g;
  Runs R;
  HD void operator()(i64 r) const {
    const i32 m = g.rowmask[r], y = i32(r - g.rs[m]), W = g.W[m], H = g.H[m];
    const i32 sent = R.rb[g.rs[m]];
    const i32 b_lo = R.rb[r] + (y == 0), b_hi = R.rb[r + 1];
    for (i32 b = b_lo; b < b_hi; ++b)
      if (y == 0 || y == H - 1 || R.bx0[b] == 0 || R.bx1[b] == W - 1) unite(R.bpar, b, sent);
    if (y == 0) return;
#ifdef K1_MUTATE_FG4  // test build (D-G16): foreground 4-adjacent, background 8-adjacent
    constexpr i32 kFg = 0, kBg = 1;
#else
    constexpr i32 kFg = 1, kBg = 0;  // foreground 8-adjacent: a.x0 <= b.x1 + 1; background 4-adjacent
#endif
    for (i32 a = R.rf[r], ae = R.rf[r + 1], b = R.rf[r - 1], bend = R.rf[r]; a < ae && b < bend;) {
      if (R.fx0[a] <= R.fx1[b] + kFg && R.fx0[b] <= R.fx1[a] + kFg) unite(R.fpar, a, b);
      const i32 xa = R.fx1[a], xb = R.fx1[b];
      if (xa <= xb) ++a;
      if (xb <= xa) ++b;
    }
    for (i32 a = b_lo, ae = b_hi, b = R.rb[r - 1] + (y - 1 == 0), bend = R.rb[r]; a < ae && b < bend;) {
      if (R.bx0[a] <= R.bx1[b] + kBg && R.bx0[b] <= R.bx1[a] + kBg) unite(R.bpar, a, b);
      const i32 xa = R.bx1[a], xb = R.bx1[b];
      if (xa <= xb) ++a;
      if (xb <= xa) ++b;
    }
  }
};

struct Flatten {
  Runs R;
  i32 nf;
  HD void operator()(i64 i) const {
    if (i < nf) R.fpar[i] = find(R.fpar, i32(i));
    else R.bpar[i - nf] = find(R.bpar, i32(i - nf));
  }
};

// ---------------------------------------------------------------- rings: order, starts, parts
// Order key of each hole ring: its first east edge in raster order. An east edge of the hole ring of B is the top
// side of a foreground pixel p of S(B) (the component surrounding B) whose upper neighbour is in B. S(B) is the
// component of the foreground run just left of B's first run (DESIGN.md lemma 3).
struct HoleKey {
  Geo g;
  Runs R;
  i32* err;
  HD void operator()(i64 r) const {
    const i32 m = g.rowmask[r], y = i32(r - g.rs[m]), W = g.W[m];
    if (y == 0) return;
    for (i32 a = R.rf[r], ae = R.rf[r + 1], b = R.rb[r - 1] + (y - 1 == 0), bend = R.rb[r]; a < ae && b < bend;) {
      if (R.fx0[a] <= R.bx1[b] && R.bx0[b] <= R.fx1[a]) {
        const i32 B = R.bpar[b];
        if (R.bx0[B] >= 0) {  // a hole (not the outside)
          const i32 left = R.bleft[B];
          if (left < 0) set_flag(err);
          else if (R.fpar[a] == R.fpar[left]) atomic_min(R.hkey + B, y * W + mx(R.fx0[a], R.bx0[b]));
        }
      }
      const i32 xa = R.fx1[a], xb = R.bx1[b];
      if (xa <= xb) ++a;
      if (xb <= xa) ++b;
    }
  }
};

// One bit per ring in the order bitmap (at its order pixel) and in the start bitmap (at its start pixel: a
// foreground pixel for an outer ring, a background one for a hole); the mask's Betti numbers (one per root).
struct RingMark {
  Geo g;
  Runs R;
  i32 nf;
  u32 *obits, *sbits;
  i64* betti;
  i32* err;
  HD void operator()(i64 i) const {
    if (i < nf) {
      if (R.fpar[i] != i) return;
      const i32 r = R.frow[i], m = g.rowmask[r], y = r - g.rs[m], x = R.fx0[i];
      atomic_or(obits + g.bitword(m, x, y), 1u << (x & 31));
      atomic_or(sbits + g.bitword(m, x, y), 1u << (x & 31));
      atomic_add(betti + 2 * i64(m), i64(1));
    } else {
      const i32 b = i32(i - nf);
      if (R.bpar[b] != b || R.bx0[b] < 0) return;
      const i32 r = R.brow[b], m = g.rowmask[r], y = r - g.rs[m], x = R.bx0[b], W = g.W[m];
      atomic_or(sbits + g.bitword(m, x, y), 1u << (x & 31));
      atomic_add(betti + 2 * i64(m) + 1, i64(1));
      const i32 k = R.hkey[b];
      if (k == INT_MAX) {
        set_flag(err);
        return;
      }
      const i32 ky = k / W, kx = k - ky * W;
      atomic_or(obits + g.bitword(m, kx, ky), 1u << (kx & 31));
    }
  }
};

struct Popc {  // set bits per word of `bits` (and of `also`, when given)
  const u32 *bits, *also;
  i32* cnt;
  HD void operator()(i64 w) const { cnt[w] = popc(bits[w] & (also ? also[w] : ~0u)); }
};

HD inline i32 rank_of(const u32* bits, const i32* scan, i32 w, i32 x) {
  return scan[w] + popc(bits[w] & ((1u << (x & 31)) - 1u));
}

struct Seeds {
  i32 *x, *y, *d, *m;  // start edge of each seed; seeds [0, R) are the ring starts, in ring order
  i32 *cnt, *next, *off, *ring, *dead;
  i64* area;
};

struct RingInit {
  Geo g;
  Runs R;
  i32 nf;
  const u32 *obits, *sbits, *words;
  const i32 *oscan, *sscan, *pscan;  // pscan: outer starts (start bits on foreground pixels) before each word
  i32 *ring_of_start, *part;
  Seeds S;
  HD void operator()(i64 i) const {
    i32 r, x, ox, oy, d, f;  // f: the root foreground run of the ring's part
    if (i < nf) {
      if (R.fpar[i] != i) return;
      r = R.frow[i], x = R.fx0[i], d = 0, f = i32(i);
    } else {
      const i32 b = i32(i - nf);
      if (R.bpar[b] != b || R.bx0[b] < 0) return;
      r = R.brow[b], x = R.bx0[b], d = 1, f = R.bleft[b] < 0 ? -1 : R.fpar[R.bleft[b]];
    }
    const i32 m = g.rowmask[r], y = r - g.rs[m];
    if (d == 0) ox = x, oy = y;
    else oy = R.hkey[i - nf] / g.W[m], ox = R.hkey[i - nf] - oy * g.W[m];
    const i32 ring = rank_of(obits, oscan, g.bitword(m, ox, oy), ox);
    ring_of_start[rank_of(sbits, sscan, g.bitword(m, x, y), x)] = ring;
    S.x[ring] = x, S.y[ring] = y, S.d[ring] = d, S.m[ring] = m;
    if (f < 0) {  // HoleKey has flagged it
      part[ring] = -1;
      return;
    }
    const i32 fx = R.fx0[f], fw = g.bitword(m, fx, R.frow[f] - g.rs[m]);
    part[ring] = pscan[fw] - pscan[g.ws[m]] + popc(sbits[fw] & words[fw] & ((1u << (fx & 31)) - 1u));
  }
};

struct SeedFlag {
  Geo g;
  Runs R;
  i32* flag;
  HD void operator()(i64 i) const {
    const i32 r = R.frow[i], m = g.rowmask[r];
    flag[i] = ((r - g.rs[m]) % g.K) == 0;
  }
};

struct SeedEmit {
  Geo g;
  Runs R;
  const i32 *flag, *soff, *counts;  // counts[0] = rings
  i32* seedmap;
  Seeds S;
  HD void operator()(i64 i) const {
    if (!flag[i]) return;
    const i32 r = R.frow[i], m = g.rowmask[r], y = r - g.rs[m], W = g.W[m];
    const i32 id = counts[0] + 2 * soff[i], x0 = R.fx0[i], x1 = R.fx1[i];
    i32* row = seedmap + g.ss[m] + (y / g.K) * (W + 1);
    row[x0] = id, row[x1 + 1] = id + 1;
    S.x[id] = x0, S.y[id] = y + 1, S.d[id] = 3, S.m[id] = m;                  // left side of pixel (x0, y), north
    S.x[id + 1] = x1 + 1, S.y[id + 1] = y, S.d[id + 1] = 1, S.m[id + 1] = m;  // right side of pixel (x1, y), south
  }
};

struct SetCounts {  // counts[0] = rings, [1] = seeds; Header sets [2] = corners; [4] = the error flag
  const i32 *oscan, *soff;
  i32 nw, nf, K;
  i32* counts;
  HD void operator()(i64) const {
    counts[0] = oscan[nw];
    counts[1] = counts[0] + (K > 0 ? 2 * soff[nf] : 0);
  }
};

// ---------------------------------------------------------------- tracing
// Walk boundary edges from a seed to the next seed. pass 1: corners, area, next seed; pass 2: write the corners.
struct Trace {
  Geo g;
  const u32 *words, *sbits;
  const i32 *sscan, *ring_of_start, *seedmap, *counts, *rbase;
  Seeds S;
  i32* pts;
  i32* err;
  int pass;

  HD i32 seed_of(i32 m, i32 x, i32 y, int d) const {  // seed id of edge (x, y, d), or -1
    if (d <= 1 && x < g.W[m] && y < g.H[m]) {
      const i32 w = g.bitword(m, x, y);
      if ((sbits[w] >> (x & 31)) & 1u) return ring_of_start[rank_of(sbits, sscan, w, x)];
    }
    if (g.K > 0) {
      if (d == 1 && y % g.K == 0) return seedmap[g.ss[m] + (y / g.K) * (g.W[m] + 1) + x];
      if (d == 3 && (y - 1) % g.K == 0) return seedmap[g.ss[m] + ((y - 1) / g.K) * (g.W[m] + 1) + x];
    }
    return -1;
  }

  HD void operator()(i64 i) const {
    const i32 s = i32(i);
    if (s >= counts[1]) return;
    if (pass == 2 && S.dead[s]) return;
    const i32 m = S.m[s];
    i32 x = S.x[s], y = S.y[s];
    int d = S.d[s];
    if (pass == 1 && s >= counts[0] && seed_of(m, x, y, d) != s) {  // a row seed on a ring start: the start covers it
      S.dead[s] = 1;
      return;
    }
    int prev = g.pred_dir(words, m, x, y, d);
    i64 k = pass == 2 ? i64(rbase[S.ring[s]]) + S.off[s] : 0;
    i32 corners = 0, nxt = -1;
    i64 area = 0;
    const i64 cap = 4 * (i64(g.W[m]) + 1) * (i64(g.H[m]) + 1);
    for (i64 step = 0; nxt < 0; ++step) {
      if (step > cap || d < 0) {
        set_flag(err);
        return;
      }
      if (d != prev) {
        if (pass == 2) pts[2 * k] = x, pts[2 * k + 1] = y, ++k;
        ++corners;
      }
      area += area_term(x, y, d);
      x += step_x(d), y += step_y(d);
      prev = d;
      d = g.succ_dir(words, m, x, y, d);
      if (d < 0) {
        set_flag(err);
        return;
      }
      nxt = seed_of(m, x, y, d);
    }
    if (pass == 1) S.cnt[s] = corners, S.area[s] = area, S.next[s] = nxt;
  }
};

struct RingWalk {  // per ring: follow its seeds; each segment gets the offset of its first corner
  const i32* counts;
  Seeds S;
  i32* rcount;
  i64* rarea;
  i32* err;
  i32 cap;
  HD void operator()(i64 i) const {
    const i32 g0 = i32(i);
    if (g0 >= counts[0]) return;
    i32 s = g0, off = 0, n = 0;
    i64 a = 0;
    do {
      S.off[s] = off, S.ring[s] = g0;
      off += S.cnt[s], a += S.area[s];
      s = S.next[s];
      if (++n > cap || s < 0) {
        set_flag(err);
        return;
      }
    } while (s != g0);
    rcount[g0] = off, rarea[g0] = a;
  }
};

struct Header {  // counts[2] = corners in all; the first ring of each mask (the order scan at its first word)
  Geo g;
  const i32 *oscan, *rbase;
  i32 rmax;
  i32 *counts, *ring_start;
  HD void operator()(i64 m) const {
    if (m == 0) counts[2] = rbase[rmax];
    ring_start[m] = oscan[g.ws[m]];
  }
};

// ---------------------------------------------------------------- one chunk
// Masks [m0, m1) on one lane: layout A (inputs, words, bitmaps) in L.dev, layout B (runs, rings, seeds, corners) in
// L.score, CUB's scratch in L.scan, three host waits. Returns false (nothing appended) when layout B would not fit
// (above max_scratch_bytes, or past the i32 counts) and the chunk has more than one mask, so the caller halves it
// (the steps before wait 1 are cheap to redo); one mask that does not fit raises.
bool mask_chunk(const Masks& in, i64 m0, i64 m1, const MaskOptions& opt, MaskResult& out, Lane& L, Timer& t) {
  const i32 n = i32(m1 - m0), K = opt.seed_rows;
  std::vector<i32> W(n), H(n), wpr(n), ws(n + 1, 0), rs(n + 1, 0), vs(n + 1, 0), ss(n + 1, 0);
  std::vector<i64> start(n + 1, 0);
  for (i32 k = 0; k < n; ++k) {  // i32 counts: a chunk holds at most kMaxChunkPixels + one mask (check_masks)
    W[k] = in.W[m0 + k], H[k] = in.H[m0 + k], wpr[k] = (W[k] + 31) / 32;
    ws[k + 1] = ws[k] + H[k] * wpr[k];
    rs[k + 1] = rs[k] + H[k];
    vs[k + 1] = vs[k] + H[k] + 1;
    ss[k + 1] = ss[k] + (K > 0 ? (H[k] + K - 1) / K * (W[k] + 1) : 0);
    start[k + 1] = in.start[m0 + k + 1] - in.start[m0];
  }
  const i64 npx = start[n];
  const i32 NW = ws[n], NR = rs[n], NV = vs[n];
  std::vector<i32> rowmask(NR);
  for (i32 k = 0; k < n; ++k) std::fill(rowmask.begin() + rs[k], rowmask.begin() + rs[k + 1], k);

  Layout la;  // ---- layout A
  const size_t o_bytes = la.put<u8>(npx), o_start = la.put<i64>(n + 1), o_W = la.put<i32>(n), o_H = la.put<i32>(n),
               o_wpr = la.put<i32>(n), o_ws = la.put<i32>(n + 1), o_rs = la.put<i32>(n + 1), o_vs = la.put<i32>(n + 1),
               o_ss = la.put<i32>(n + 1), o_rowmask = la.put<i32>(NR);
  const size_t in_bytes = la.size;
  const size_t o_words = la.put<u32>(NW), o_cf = la.put<i32>(NR + 1), o_cb = la.put<i32>(NR + 1),
               o_rf = la.put<i32>(NR + 1), o_rb = la.put<i32>(NR + 1), o_ocnt = la.put<i32>(NW + 1),
               o_scnt = la.put<i32>(NW + 1), o_pcnt = la.put<i32>(NW + 1), o_oscan = la.put<i32>(NW + 1),
               o_sscan = la.put<i32>(NW + 1), o_pscan = la.put<i32>(NW + 1), o_seedmap = la.put<i32>(ss[n]);
  const size_t z0 = la.size;  // [z0, end) zeroed
  const size_t o_obits = la.put<u32>(NW), o_sbits = la.put<u32>(NW), o_quads = la.put<i64>(5 * i64(n)),
               o_betti = la.put<i64>(2 * i64(n)), o_counts = la.put<i32>(8);
  u8* A = L.dev.get(la.size);
  u8* h = L.in.get(in_bytes);
  std::copy_n(in.bits + in.start[m0], npx, ptr<u8>(h, o_bytes));  // copy_n: a count of 0 never reads
  std::copy_n(start.data(), n + 1, ptr<i64>(h, o_start));
  std::copy_n(W.data(), n, ptr<i32>(h, o_W));
  std::copy_n(H.data(), n, ptr<i32>(h, o_H));
  std::copy_n(wpr.data(), n, ptr<i32>(h, o_wpr));
  std::copy_n(ws.data(), n + 1, ptr<i32>(h, o_ws));
  std::copy_n(rs.data(), n + 1, ptr<i32>(h, o_rs));
  std::copy_n(vs.data(), n + 1, ptr<i32>(h, o_vs));
  std::copy_n(ss.data(), n + 1, ptr<i32>(h, o_ss));
  std::copy_n(rowmask.data(), NR, ptr<i32>(h, o_rowmask));

  const Stream st = L.stream;
  const Geo g{ptr<i32>(A, o_W), ptr<i32>(A, o_H), ptr<i32>(A, o_wpr), ptr<i32>(A, o_ws), ptr<i32>(A, o_rs),
              ptr<i32>(A, o_vs), ptr<i32>(A, o_ss), ptr<i32>(A, o_rowmask), n, K};
  u32 *words = ptr<u32>(A, o_words), *obits = ptr<u32>(A, o_obits), *sbits = ptr<u32>(A, o_sbits);
  i32 *counts = ptr<i32>(A, o_counts), *err = counts + 4;
  i32 *oscan = ptr<i32>(A, o_oscan), *sscan = ptr<i32>(A, o_sscan), *pscan = ptr<i32>(A, o_pscan);
#ifdef GC_ALIAS_SCAN_FOR_TEST
  Block<true>& scan1 = L.dev;  // planted defect (tests/variants.py "alias"): CUB scratch on layout A
#else
  Block<true>& scan1 = L.scan;
#endif
  u8* o = L.out.get(64);
  t.stage("ccl", st, [&] {
    to_device(A, h, in_bytes, st);
    zero(A + z0, la.size - z0, st);
    for (size_t c : {o_cf, o_cb}) zero(A + c, sizeof(i32) * size_t(NR + 1), st);
    for (size_t c : {o_ocnt, o_scnt, o_pcnt}) zero(A + c, sizeof(i32) * size_t(NW + 1), st);
    pack_words(st, Pack{g, ptr<u8>(A, o_bytes), ptr<i64>(A, o_start), words}, NW);
    for_each(st, NV, Quads{g, words, ptr<i64>(A, o_quads)});
    for_each(st, NR, RunCount{g, words, ptr<i32>(A, o_cf), ptr<i32>(A, o_cb)});
    scan(st, ptr<i32>(A, o_cf), ptr<i32>(A, o_rf), NR + 1, scan1);
    scan(st, ptr<i32>(A, o_cb), ptr<i32>(A, o_rb), NR + 1, L.scan);
  });
  to_host(o, A + o_rf + sizeof(i32) * NR, sizeof(i32), st);
  to_host(o + 4, A + o_rb + sizeof(i32) * NR, sizeof(i32), st);
  wait(st);  // 1: the run totals
  const i32 NF = ptr<i32>(o, 0)[0], NB = ptr<i32>(o, 4)[0];
  if (NF < 0 || NB < 0 || NF > npx || NB > npx + NR) throw std::logic_error("gpucheck: run totals out of their bound");
  const i64 rmax = i64(NF) + NB, smax = rmax + 2 * i64(NF);

  Layout lb;  // ---- layout B
  size_t o_f[4], o_b[6], o_s[8];
  for (size_t& q : o_f) q = lb.put<i32>(NF);
  for (size_t& q : o_b) q = lb.put<i32>(NB);
  const size_t o_flag = lb.put<i32>(NF + 1), o_soff = lb.put<i32>(NF + 1), o_ros = lb.put<i32>(rmax),
               o_part = lb.put<i32>(rmax), o_rbase = lb.put<i32>(rmax + 1), o_rarea = lb.put<i64>(rmax),
               o_pts = lb.put<i32>(8 * i64(NF)), o_rstart = lb.put<i32>(n + 1);
  for (size_t& q : o_s) q = lb.put<i32>(smax);
  const size_t o_sarea = lb.put<i64>(smax);
  const size_t zb = lb.size;  // [zb, end) zeroed
  const size_t o_rcount = lb.put<i32>(rmax + 1), o_dead = lb.put<i32>(smax);
  if (smax > INT_MAX / 2 || 8 * i64(NF) > INT_MAX || i64(lb.size) > opt.max_scratch_bytes) {
    if (n > 1) return false;
    throw std::invalid_argument("check_masks: one mask needs " + std::to_string(lb.size) +
                                " bytes of run scratch, above max_scratch_bytes = " +
                                std::to_string(opt.max_scratch_bytes));
  }
  u8* B = L.score.get(lb.size);
  const i32 RMAX = i32(rmax), SMAX = i32(smax);
  const Runs R{ptr<i32>(A, o_rf), ptr<i32>(A, o_rb), ptr<i32>(B, o_f[0]), ptr<i32>(B, o_f[1]), ptr<i32>(B, o_f[2]),
               ptr<i32>(B, o_f[3]), ptr<i32>(B, o_b[0]), ptr<i32>(B, o_b[1]), ptr<i32>(B, o_b[2]), ptr<i32>(B, o_b[3]),
               ptr<i32>(B, o_b[4]), ptr<i32>(B, o_b[5])};
  const Seeds S{ptr<i32>(B, o_s[0]), ptr<i32>(B, o_s[1]), ptr<i32>(B, o_s[2]), ptr<i32>(B, o_s[3]),
                ptr<i32>(B, o_s[4]), ptr<i32>(B, o_s[5]), ptr<i32>(B, o_s[6]), ptr<i32>(B, o_s[7]),
                ptr<i32>(B, o_dead), ptr<i64>(B, o_sarea)};
  i32 *ros = ptr<i32>(B, o_ros), *rcount = ptr<i32>(B, o_rcount), *rbase = ptr<i32>(B, o_rbase);
  i32 *flag = ptr<i32>(B, o_flag), *soff = ptr<i32>(B, o_soff);
  t.stage("ccl", st, [&] {
    zero(B + zb, lb.size - zb, st);
    zero(flag, sizeof(i32) * size_t(NF + 1), st);
    for_each(st, NR, RunEmit{g, words, R});
    for_each(st, NR, Merge{g, R});
    for_each(st, rmax, Flatten{R, NF});
  });
  t.stage("rings", st, [&] {
    for_each(st, NR, HoleKey{g, R, err});
    for_each(st, rmax, RingMark{g, R, NF, obits, sbits, ptr<i64>(A, o_betti), err});
    for_each(st, NW, Popc{obits, nullptr, ptr<i32>(A, o_ocnt)});
    for_each(st, NW, Popc{sbits, nullptr, ptr<i32>(A, o_scnt)});
    for_each(st, NW, Popc{sbits, words, ptr<i32>(A, o_pcnt)});
    scan(st, ptr<i32>(A, o_ocnt), oscan, NW + 1, L.scan);
    scan(st, ptr<i32>(A, o_scnt), sscan, NW + 1, L.scan);
    scan(st, ptr<i32>(A, o_pcnt), pscan, NW + 1, L.scan);
    for_each(st, rmax, RingInit{g, R, NF, obits, sbits, words, oscan, sscan, pscan, ros, ptr<i32>(B, o_part), S});
    if (K > 0) {
      for_each(st, NF, SeedFlag{g, R, flag});
      scan(st, flag, soff, NF + 1, L.scan);
    }
    for_each(st, 1, SetCounts{oscan, soff, NW, NF, K, counts});
    if (K > 0) for_each(st, NF, SeedEmit{g, R, flag, soff, counts, ptr<i32>(A, o_seedmap), S});
  });
  t.stage("trace", st, [&] {
    Trace tr{g, words, sbits, sscan, ros, ptr<i32>(A, o_seedmap), counts, rbase, S, ptr<i32>(B, o_pts), err, 1};
    for_each(st, SMAX, tr);
    for_each(st, RMAX, RingWalk{counts, S, rcount, ptr<i64>(B, o_rarea), err, SMAX});
    scan(st, rcount, rbase, rmax + 1, L.scan);
    tr.pass = 2;
    for_each(st, SMAX, tr);
    for_each(st, n + 1, Header{g, oscan, rbase, RMAX, counts, ptr<i32>(B, o_rstart)});
  });
  to_host(o, counts, sizeof(i32) * 8, st);
  wait(st);  // 2: the ring and corner totals
  i32 hdr[8];
  std::copy_n(ptr<i32>(o, 0), 8, hdr);
  if (hdr[4]) throw std::logic_error("gpucheck: polygonizer invariant violated on the device (trace, hole key, walk)");
  const i32 nr = hdr[0], ncorner = hdr[2];
  if (nr < 0 || nr > RMAX || ncorner < 0 || ncorner > 4 * NF) throw std::logic_error("gpucheck: counts out of bound");
  Layout lo;
  const size_t q_rbase = lo.put<i32>(nr + 1), q_area = lo.put<i64>(nr), q_part = lo.put<i32>(nr),
               q_pts = lo.put<i32>(2 * i64(ncorner)), q_rstart = lo.put<i32>(n + 1), q_quads = lo.put<i64>(5 * i64(n)),
               q_betti = lo.put<i64>(2 * i64(n));
  o = L.out.get(lo.size);
  to_host(o + q_rbase, rbase, sizeof(i32) * size_t(nr + 1), st);
  to_host(o + q_area, ptr<i64>(B, o_rarea), sizeof(i64) * size_t(nr), st);
  to_host(o + q_part, ptr<i32>(B, o_part), sizeof(i32) * size_t(nr), st);
  to_host(o + q_pts, ptr<i32>(B, o_pts), sizeof(i32) * 2 * size_t(ncorner), st);
  to_host(o + q_rstart, ptr<i32>(B, o_rstart), sizeof(i32) * size_t(n + 1), st);
  to_host(o + q_quads, ptr<i64>(A, o_quads), sizeof(i64) * 5 * size_t(n), st);
  to_host(o + q_betti, ptr<i64>(A, o_betti), sizeof(i64) * 2 * size_t(n), st);
  wait(st);  // 3: the outputs
  const i64 r0 = i64(out.ring_area.size()), p0 = out.pts_start.back();
  const i32 *rst = ptr<i32>(o, q_rstart), *rb = ptr<i32>(o, q_rbase), *pp = ptr<i32>(o, q_pts);
  const i64 *ar = ptr<i64>(o, q_area), *qq = ptr<i64>(o, q_quads), *bb = ptr<i64>(o, q_betti);
  for (i32 k = 1; k <= n; ++k) out.ring_start.push_back(r0 + rst[k]);
  for (i32 k = 0; k < nr; ++k) out.pts_start.push_back(p0 + rb[k + 1]), out.ring_area.push_back(ar[k] / 2);
  out.ring_part.insert(out.ring_part.end(), ptr<i32>(o, q_part), ptr<i32>(o, q_part) + nr);
  out.pts.insert(out.pts.end(), pp, pp + 2 * size_t(ncorner));
  out.quads.insert(out.quads.end(), qq, qq + 5 * size_t(n));
  out.betti.insert(out.betti.end(), bb, bb + 2 * size_t(n));
  return true;
}

void mask_range(const Masks& in, i64 m0, i64 m1, const MaskOptions& opt, MaskResult& out, Lane& L, Timer& t) {
  if (mask_chunk(in, m0, m1, opt, out, L, t)) return;
  const i64 mid = m0 + (m1 - m0) / 2;
  mask_range(in, m0, mid, opt, out, L, t);
  mask_range(in, mid, m1, opt, out, L, t);
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

// The copied-in arrays of a chunk of n shapes, nr rings and nv vertices, at the front of its layout.
void put_inputs(Layout& l, PolyChunk& c) {
  c.xy = l.put<i32>(2 * c.nv), c.roles = l.put<i32>(c.nr), c.W = l.put<i32>(c.n), c.H = l.put<i32>(c.n);
  c.ring_start = l.put<i64>(c.nr + 1), c.shape_start = l.put<i64>(c.n + 1), c.row_start = l.put<i64>(c.n + 1);
  c.word_start = l.put<i64>(c.n + 1);
  c.in_bytes = l.size;
}

PolyChunk plan(const Polygons& in, const std::vector<i64>& word_start, i64 a, i64 b) {
  PolyChunk c{};
  c.a = a, c.b = b, c.n = b - a;
  c.r0 = in.shape_start[a], c.nr = in.shape_start[b] - c.r0;
  c.v0 = in.ring_start[c.r0], c.nv = in.ring_start[c.r0 + c.nr] - c.v0;
  c.rows = 0;
  for (i64 s = a; s < b; ++s) c.rows += in.H[s];
  c.nwords = word_start[b] - word_start[a];
  Layout l;
  put_inputs(l, c);
  c.raster_wide = l.put<u8>(c.n), c.validate_wide = l.put<u8>(c.n), c.fewer = l.put<u8>(c.nr);
  c.cxy = l.put<i32>(2 * c.nv), c.len = l.put<i32>(c.nr), c.box = l.put<i32>(4 * c.nr), c.area2 = l.put<i64>(c.nr);
  c.counts = l.put<i32>(c.nr * kNumDefects), c.pairs = l.put<i32>(c.n * 2 * kPairs);  // zeroed, with defects
  c.defects = l.put<i64>(c.n * kNumDefects), c.words = l.put<u32>(c.nwords);
  c.size = l.size;
  return c;
}

// Stage bounds: one score_polygons call of 256 answers of at most 2 objects (geometry.text SCHEMAS: RIGA's disc and
// cup) at the grammar's caps (10 rings of 1,000 vertices per object). In: the chunk's copied-in arrays, about 41 MB
// (the scorer's own staging, planes and tables per object, is about 100 KB at this n). Out: the largest copy back of
// such a call: the signatures at 64 components per object (twice the largest count of the D1 replay set is 4) are
// 128 KiB, above the defects (52 KiB), the boxes (20 KiB) and the per-object results (32 KiB).
constexpr i64 kStageObjects = 256 * 2, kStageComps = 64;
size_t stage_in_bound() {
  PolyChunk c{};
  c.n = kStageObjects, c.nr = c.n * kMaxRings, c.nv = c.nr * kMaxRingVertices;
  Layout l;
  put_inputs(l, c);
  return l.size;
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
Shapes shapes(const PolyChunk& c, u8* d, i64 G) {  // the chunk's shapes in the arena at d
  return Shapes{ptr<i32>(d, c.xy),         ptr<i64>(d, c.ring_start), ptr<i64>(d, c.shape_start),
                ptr<i64>(d, c.row_start),  ptr<i64>(d, c.word_start), ptr<i32>(d, c.roles),
                ptr<i32>(d, c.W),          ptr<i32>(d, c.H),          c.n,
                G};
}

void launch(const PolyChunk& c, i64 G, Lane& L, bool masks, Timer& t) {
  u8* d = L.dev.p;
  const Stream st = L.stream;
  const Shapes s = shapes(c, d, G);
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
    for_each_warp(st, c.n, ShapeTotalsWarp{s, r, validate_wide, ptr<i64>(d, c.defects)});
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

struct Registry {  // resident tables; handle = index, never reused (a stale handle raises)
  std::vector<u8*> block;
  std::vector<size_t> bytes;
  std::vector<TableView> view;
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

// The masks bit-packed into lane L's checker arena (mask k at word word_start[k]): one upload, one wait.
const u32* upload_masks(Lane& L, const Masks& in, const std::vector<i64>& word_start) {
  const size_t bytes = sizeof(u32) * size_t(word_start.back());
  u8* h = L.in.get(bytes);
  std::fill_n(h, bytes, u8(0));
  for (i64 k = 0; k < in.n; ++k)
    pack_plane(in.bits + in.start[k], in.H[k], in.W[k], ptr<u32>(h, 0) + word_start[size_t(k)]);
  u8* d = L.dev.get(bytes);
  to_device(d, h, bytes, L.stream);
  wait(L.stream);
  return ptr<u32>(d, 0);
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

// ---------------------------------------------------------------- the pixel label chain (scorer, certify_polygons)
// A bit-plane mask of one call: H rows of words_per_row(W) words from word `words` of the call's word buffer, size-
// filtered at min_area. The box step fills its foreground box (bx, by, bw, bh), its area and its crop's offset.
struct Plane {
  i64 words, min_area, area, crop;
  i32 W, H, bx, by, bw, bh;
};

// Per frame row: first and last foreground column, foreground pixels.
struct RowBox {
  const Plane* pl;
  const i64* row_start;
  i64 n;
  const u32* words;
  i32* rbox;
  HD void operator()(i64 i) const {
    const i64 o = segment(row_start, n, i);
    const i32 wpr = words_per_row(pl[o].W);
    const u32* w = words + pl[o].words + (i - row_start[o]) * wpr;
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

// Per plane: box (x, y, w, h) and area.
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

// Union-find labeling on each plane's crop: its box padded by one pixel (8-connected foreground, 4-connected
// background; the outside is the background component of pixel 0).
struct Crops {  // the call's crops, one item per crop row; pixel p of a crop = y (bw + 2) + x, local to it
  const Plane* pl;
  const i64* crow_start;
  i64 n;
  u8* img;             // 1 = foreground
  i32 *fg, *bg, *aux;  // labels (root = the component's first raster pixel), per-root counts
  i32* ncomp;          // per plane: components after the filter
  i32* changed;        // per plane: rows the size filter changed (0 = the first labels stay)
  i64* relabelled;     // per plane: crop rows the second label pair relabelled (the label_pairs counter)
  HD i64 row(i64 r, i32& y) const {
    const i64 o = segment(crow_start, n, r);
    y = i32(r - crow_start[o]);
    return o;
  }
};

struct CropFill {
  Crops c;
  const u32* words;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Plane& b = c.pl[c.row(r, y)];
    const i32 cw = b.bw + 2, wpr = words_per_row(b.W);
    u8* img = c.img + b.crop + i64(y) * cw;
    const bool inside = y >= 1 && y <= b.bh;
    const u32* row = words + b.words + i64(b.by + (inside ? y - 1 : 0)) * wpr;
    for (i32 x = 0; x < cw; ++x) {
      const i32 fx = b.bx + x - 1;
      img[x] = (inside && x >= 1 && x <= b.bw) ? u8(row[fx / 32] >> (fx % 32) & 1) : u8(0);
    }
  }
};

// The three labeling steps of one class: which = 1 (foreground, 8-connected) or 0 (background, 4-connected). With
// `lazy` (the second pair, after the size filter) a row is labelled again only if the filter changed its plane: else
// the first labels are those of the same image (roots are first pixels, so labels do not depend on the run).
struct LabelInit {
  Crops c;
  int which;
  bool lazy;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const i64 o = c.row(r, y);
    if (lazy && !c.changed[o]) return;
    if (lazy && which) atomic_add(c.relabelled + o, 1);
    const Plane& b = c.pl[o];
    const i32 cw = b.bw + 2;
    i32* par = (which ? c.fg : c.bg) + b.crop;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p) par[p] = c.img[b.crop + p] == which ? p : -1;
  }
};

struct LabelMerge {
  Crops c;
  int which;
  bool lazy;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const i64 o = c.row(r, y);
    if (lazy && !c.changed[o]) return;
    const Plane& b = c.pl[o];
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
  bool lazy;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const i64 o = c.row(r, y);
    if (lazy && !c.changed[o]) return;
    const Plane& b = c.pl[o];
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
    const Plane& b = c.pl[c.row(r, y)];
    const i32 cw = b.bw + 2;
    const i32* par = (which ? c.fg : c.bg) + b.crop;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p)
      if (c.img[b.crop + p] == which) atomic_add(c.aux + b.crop + par[p], 1);
  }
};

// Size filter (topology.size_filter): which = 1 drops foreground components below min_area, which = 0 fills holes
// (background components other than the outside, root 0) below min_area. aux holds the class's counts. A row with a
// changed pixel adds 1 to its plane's `changed`.
struct Filter {
  Crops c;
  int which;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const i64 o = c.row(r, y);
    const Plane& b = c.pl[o];
    const i32 cw = b.bw + 2;
    const i32* par = (which ? c.fg : c.bg) + b.crop;
    bool hit = false;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p) {
      if (c.img[b.crop + p] != which || (!which && par[p] == 0) || c.aux[b.crop + par[p]] >= b.min_area) continue;
      c.img[b.crop + p] = u8(1 - which), hit = true;
    }
    if (hit) atomic_add(c.changed + o, 1);
  }
};

// A hole's first raster pixel is its root; the pixel above it is foreground of the enclosing component (an island
// inside the hole lies strictly below the hole's top row): to[foreground root] += 1 (to: per crop pixel, zeroed).
struct Holes {
  Crops c;
  i32* to;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Plane& b = c.pl[c.row(r, y)];
    const i32 cw = b.bw + 2;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p)
      if (p != 0 && c.img[b.crop + p] == 0 && c.bg[b.crop + p] == p)
        atomic_add(to + b.crop + c.fg[b.crop + p - cw], 1);
  }
};

// Every component's hole count into the plane's list: the first ncomp entries of its crop in bg (no longer needed).
struct Gather {
  Crops c;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const i64 o = c.row(r, y);
    const Plane& b = c.pl[o];
    const i32 cw = b.bw + 2;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p)
      if (c.img[b.crop + p] == 1 && c.fg[b.crop + p] == p) c.bg[b.crop + atomic_inc(c.ncomp + o)] = c.aux[b.crop + p];
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

struct SortList {  // per plane: its list, descending (the signature)
  Crops c;
  HD void operator()(i64 o) const { sort_descending(c.bg + c.pl[o].crop, c.ncomp[o]); }
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

struct SigCompact {  // the lists, one after the other, into aux (no longer needed)
  Crops c;
  const i64* sig_start;
  HD void operator()(i64 o) const {
    for (i32 j = 0; j < c.ncomp[o]; ++j) c.aux[sig_start[o] + j] = c.bg[c.pl[o].crop + j];
  }
};

// The label chain's arrays in a call's layout: the planes and their crop rows (copied in, at the front), then per
// crop pixel the image, the labels and the counts.
struct CropLayout {
  i64 n, crop;
  size_t pl = 0, crows = 0, img = 0, fg = 0, bg = 0, aux = 0, nc = 0, ch = 0, rl = 0, end = 0;
  void front(Layout& l) { pl = l.put<Plane>(n), crows = l.put<i64>(n + 1); }
  void back(Layout& l) {  // the per-plane counters last and together: label_chain zeroes [nc, end) at once
    img = l.put<u8>(crop), fg = l.put<i32>(crop), bg = l.put<i32>(crop), aux = l.put<i32>(crop);
    nc = l.put<i32>(n), ch = l.put<i32>(n), rl = l.put<i64>(n), end = l.size;
  }
  void stage(u8* h, const std::vector<Plane>& p, const std::vector<i64>& crow_start) const {
    std::copy_n(p.data(), n, ptr<Plane>(h, pl));
    std::copy_n(crow_start.data(), n + 1, ptr<i64>(h, crows));
  }
  Crops at(u8* d) const {
    return Crops{ptr<Plane>(d, pl), ptr<i64>(d, crows), n,        ptr<u8>(d, img),
                 ptr<i32>(d, fg),   ptr<i32>(d, bg),    ptr<i32>(d, aux), ptr<i32>(d, nc),
                 ptr<i32>(d, ch),   ptr<i64>(d, rl)};
  }
};

// Round trip 1 of a label-chain call: every plane's foreground box and area (one host wait), then its crop offset
// and crop rows. Returns the crop layout of the call (crop rows in crow_start). Uses the head of the scorer arena.
CropLayout boxes(Lane& L, std::vector<Plane>& pl, const u32* words, std::vector<i64>& crow_start, Timer& t) {
  const i64 n = i64(pl.size());
  const Stream st = L.stream;
  std::vector<i64> row_start(size_t(n + 1), 0);
  for (i64 o = 0; o < n; ++o) row_start[size_t(o + 1)] = row_start[size_t(o)] + pl[size_t(o)].H;
  const i64 rows = row_start[size_t(n)];
  Layout a;
  const size_t a_pl = a.put<Plane>(n), a_rows = a.put<i64>(n + 1), a_rbox = a.put<i32>(3 * rows),
               a_box = a.put<i64>(5 * n);
  u8* d = L.score.get(a.size);
  u8* h = L.in.get(a_rbox);
  std::copy_n(pl.data(), n, ptr<Plane>(h, a_pl));
  std::copy_n(row_start.data(), n + 1, ptr<i64>(h, a_rows));
  to_device(d, h, a_rbox, st);
  t.stage("label", st, [&] {  // the crop boxes, which the labeling works on
    for_each(st, rows, RowBox{ptr<Plane>(d, a_pl), ptr<i64>(d, a_rows), n, words, ptr<i32>(d, a_rbox)});
    for_each(st, n, ObjBox{ptr<i64>(d, a_rows), ptr<i32>(d, a_rbox), ptr<i64>(d, a_box)});
  });
  const i64* box = ptr<i64>(L.out.get(sizeof(i64) * size_t(5 * n)), 0);
  to_host(L.out.p, d + a_box, sizeof(i64) * size_t(5 * n), st);
  wait(st);
  crow_start.assign(size_t(n + 1), 0);
  i64 crop = 0;
  for (i64 o = 0; o < n; ++o) {
    Plane& b = pl[size_t(o)];
    const i64* x = box + 5 * o;
    b.bx = i32(x[0]), b.by = i32(x[1]), b.bw = i32(x[2]), b.bh = i32(x[3]), b.area = x[4], b.crop = crop;
    if (b.area) crop += i64(b.bw + 2) * (b.bh + 2);
    crow_start[size_t(o + 1)] = crow_start[size_t(o)] + (b.area ? b.bh + 2 : 0);
  }
  return CropLayout{n, crop};
}

// Queues the chain on the crops: fill them from the words, label, size-filter the foreground then the background
// (topology.size_filter), label again the planes the filter changed (topology.py's skip). Stages "label", "count" and
// "filter". `counters` = [ncomp, end of relabelled) of the layout: zeroed first, as the arenas are never cleared.
void label_chain(Stream st, const Crops& c, i64 crows, i64 crop, const u32* words, size_t counters, Timer& t) {
  auto label = [&](int which, bool lazy) {
    t.stage("label", st, [&] {
      for_each(st, crows, LabelInit{c, which, lazy});
      for_each(st, crows, LabelMerge{c, which, lazy});
      for_each(st, crows, LabelFlatten{c, which, lazy});
    });
  };
  auto count = [&](int which) {
    t.stage("count", st, [&] {
      zero(c.aux, sizeof(i32) * size_t(crop), st);
      for_each(st, crows, Count{c, which});
    });
  };
  auto filter = [&](int which) { t.stage("filter", st, [&] { for_each(st, crows, Filter{c, which}); }); };
  zero(c.ncomp, counters, st);
  t.stage("label", st, [&] { for_each(st, crows, CropFill{c, words}); });
  label(1, false), count(1), filter(1);
  label(0, false), count(0), filter(0);
  label(1, true), label(0, true);
}

// After label_chain, each plane's signature: hole counts per component, sorted descending, at the head of its crop
// in bg (ncomp of them; the labels in bg are gone).
void sig_lists(Stream st, const Crops& c, i64 crows, i64 crop, Timer& t) {
  t.stage("score", st, [&] {
    zero(c.aux, sizeof(i32) * size_t(crop), st);
    for_each(st, crows, Holes{c, c.aux});
    for_each(st, crows, Gather{c});
    for_each(st, c.n, SortList{c});
  });
}

// The signatures of all planes, after label_chain: sig_start into `ss` (device, n + 1) and the lists one after the
// other into aux. Nothing waits.
void compact_sigs(Stream st, const Crops& c, i64* ss, Timer& t) {
  t.stage("score", st, [&] {
    for_each(st, 1, SigScan{c.ncomp, c.n, ss});
    for_each(st, c.n, SigCompact{c, ss});
  });
}

// The last round trip of a label-chain call: the `total` signature entries that compact_sigs left in aux.
std::vector<i32> read_sigs(Lane& L, const Crops& c, i64 total) {
  std::vector<i32> sig(static_cast<size_t>(total));
  to_host(L.out.get(sizeof(i32) * size_t(total)), c.aux, sizeof(i32) * size_t(total), L.stream);
  wait(L.stream);
  std::copy_n(ptr<i32>(L.out.p, 0), total, sig.data());
  return sig;
}

// ---------------------------------------------------------------- terms
constexpr int kPart = 7;  // row partials after the k rater intersections: vote intersection, spill count and sum, miss
                          // count and sum (sums in units of 2^-23), column and row index sums of the mask's pixels

struct Obj {  // one object of a scoring call: its tables and its row partials (t.H rows of t.k + kPart)
  TableView t;
  i64 part;
};

struct ObjOut {  // copied back per object (8-byte fields only: no padding for initcheck to flag)
  i64 area, ncomp, relabelled, sum_x, sum_y;
  double iou_tc, f5, g_area, topo_freq;
};

size_t stage_out_bound() {
  const size_t n = kStageObjects, objs = (sizeof(ObjOut) * n + 255) / 256 * 256 + sizeof(i64) * (n + 1);
  return std::max({sizeof(i64) * size_t(kNumDefects) * n, sizeof(i64) * 5 * n, objs,
                   sizeof(i32) * n * size_t(kStageComps)});
}

HD inline i64 units(float d) { return i64(double(d) * 8388608.0); }  // d in units of 2^-23: exact for 0 and d >= 1

// Per frame row: intersections with each rater and the vote, spill (mask outside the union) and core miss (core outside
// the mask) pixel counts and distance sums, and the centroid sums: sx = column indices of the mask's pixels, per word w
// in closed form (bit b = column 32 w + b, the sum of b over set bits is sum_k 2^k popc(m & bit k of the index set)),
// sy = y x pixels. The bits past W are 0 in every plane and mask.
struct ScoreRows {
  const Obj* obj;
  const Plane* pl;
  const i64* row_start;
  i64 n;
  const u32* words;
  i64* parts;
  HD void operator()(i64 i) const {
    const i64 o = segment(row_start, n, i);
    const Obj& b = obj[o];
    const TableView& t = b.t;
    const i64 y = i - row_start[o], wpr = words_per_row(t.W);
    const u32* m = words + pl[o].words + y * wpr;
    i64* out = parts + b.part + y * (t.k + kPart);
    for (i32 j = 0; j < t.k; ++j) {
      const u32* r = t.raters + (j * i64(t.H) + y) * wpr;
      i64 s = 0;
      for (i64 w = 0; w < wpr; ++w) s += popc(m[w] & r[w]);
      out[j] = s;
    }
    i64 vote = 0, spill_n = 0, spill = 0, miss_n = 0, miss = 0, sx = 0, px = 0;
    const float* d = t.dist + y * t.W;
    for (i64 w = 0; w < wpr; ++w) {
      const i64 n = popc(m[w]);
      sx += 32 * w * n + popc(m[w] & 0xAAAAAAAAu) + 2 * popc(m[w] & 0xCCCCCCCCu) + 4 * popc(m[w] & 0xF0F0F0F0u) +
            8 * popc(m[w] & 0xFF00FF00u) + 16 * popc(m[w] & 0xFFFF0000u);
      px += n;
      vote += popc(m[w] & t.vote[y * wpr + w]);
      u32 s = m[w] & ~t.uni[y * wpr + w], c = t.core[y * wpr + w] & ~m[w];
      spill_n += popc(s), miss_n += popc(c);
      for (; s; s &= s - 1) spill += units(d[32 * w + ctz(s)]);
      for (; c; c &= c - 1) miss += units(d[32 * w + ctz(c)]);
    }
    i64* f = out + t.k;
    f[0] = vote, f[1] = spill_n, f[2] = spill, f[3] = miss_n, f[4] = miss, f[5] = sx, f[6] = y * px;
  }
};

HD inline double iou(i64 inter, i64 a, i64 b) {
  const i64 uni = a + b - inter;
  return uni ? double(inter) / double(uni) : 1.0;
}

HD inline double sum_value(i128 units) {  // exact below 2^53 units (2^30 px of distance), deterministic above
  const i128 lim = i128(1) << 62;
  if (units < lim && units > -lim) return double(i64(units)) / 8388608.0;
  return (double(i64(units >> 32)) * 4294967296.0 + double(u32(units & 0xFFFFFFFF))) / 8388608.0;
}

// Per object: the terms from its sorted signature and row partials (score_mask's formulas, operation by operation).
struct ObjTerms {
  const Obj* obj;
  Crops c;
  const i64* parts;
  ObjOut* out;
  HD static bool same_sig(const TableView& t, i32 r, const i32* s, i32 m) {
    if (t.sig_start[r + 1] - t.sig_start[r] != m) return false;
    for (i32 j = 0; j < m; ++j)
      if (t.sig[t.sig_start[r] + j] != s[j]) return false;
    return true;
  }
  HD void operator()(i64 o) const {
    const Obj& b = obj[o];
    const Plane& p = c.pl[o];
    const TableView& t = b.t;
    const i32 m = c.ncomp[o], F = t.k + kPart;
    const i32* s = c.bg + p.crop;
    const i64* rows = parts + b.part;
    i64 vote = 0, spill_n = 0, miss_n = 0, sx = 0, sy = 0;
    i128 spill = 0, miss = 0;
    for (i32 y = 0; y < t.H; ++y) {
      const i64* f = rows + i64(y) * F + t.k;
      vote += f[0], spill_n += f[1], spill += f[2], miss_n += f[3], miss += f[4], sx += f[5], sy += f[6];
    }
    i32 same = 0, voters = 0;
    double best = 0;
    for (i32 r = 0; r < t.k; ++r) {
      if (!t.voter[r]) continue;
      ++voters;
      if (!same_sig(t, r, s, m)) continue;
      i64 inter = 0;
      for (i32 y = 0; y < t.H; ++y) inter += rows[i64(y) * F + r];
      const double v = iou(inter, p.area, t.rater_area[r]);
      if (!same || v > best) best = v;  // Python's max: the first of equal values
      ++same;
    }
    const double d = (sum_value(spill) - 0.5 * double(spill_n)) + (sum_value(miss) - 0.5 * double(miss_n));
    double g = 0.0;
    if (p.area > 0) {
      i64 lo = t.rater_area[0], hi = t.rater_area[0];
      for (i32 r = 1; r < t.k; ++r) lo = mn(lo, t.rater_area[r]), hi = mx(hi, t.rater_area[r]);
      const double up = double(hi) / double(p.area), down = lo > 0 ? double(p.area) / double(lo) : 1.0;
      g = (up < 1.0 ? up : 1.0) * (down < 1.0 ? down : 1.0);
    }
    out[o] = ObjOut{p.area, m, c.relabelled[o], sx, sy, same ? best : iou(vote, p.area, t.n_vote),
                    1.0 / (1.0 + (2.0 * d / double(t.perimeter)) / (t.t5 * t.t5)), g, double(same) / double(voters)};
  }
};

// Scores the masks at `words` in lane L's checker arena (mask s at word_start[s]) against their tables. The staging
// buffers and the scorer's arena are free (the caller waited for the lane).
void score_words(Lane& L, const std::vector<TableView>& tv, const u32* words, const std::vector<i64>& word_start,
                 Scores& out, Timer& t) {
  const i64 n = i64(tv.size());
  const Stream st = L.stream;
  std::vector<Plane> pl(static_cast<size_t>(n));
  std::vector<Obj> obj(static_cast<size_t>(n));
  std::vector<i64> row_start(size_t(n + 1), 0), crow_start;
  i64 parts = 0;
  for (i64 o = 0; o < n; ++o) {
    const TableView& v = tv[size_t(o)];
    pl[size_t(o)] = Plane{word_start[size_t(o)], v.min_area, 0, 0, v.W, v.H, 0, 0, 0, 0};
    obj[size_t(o)] = Obj{v, parts};
    parts += i64(v.H) * (v.k + kPart);
    row_start[size_t(o + 1)] = row_start[size_t(o)] + v.H;
  }
  const i64 rows = row_start[size_t(n)];
  CropLayout cl = boxes(L, pl, words, crow_start, t);
  const i64 crows = crow_start[size_t(n)];

  // Labeling, row partials, terms.
  Layout l;
  cl.front(l);
  const size_t o_obj = l.put<Obj>(n), o_rows = l.put<i64>(n + 1), up = l.size, o_parts = l.put<i64>(parts);
  cl.back(l);
  const size_t o_out = l.put<ObjOut>(n), o_ss = l.put<i64>(n + 1);
  u8* d = L.score.get(l.size);
  u8* h = L.in.get(up);
  cl.stage(h, pl, crow_start);
  std::copy_n(obj.data(), n, ptr<Obj>(h, o_obj));
  std::copy_n(row_start.data(), n + 1, ptr<i64>(h, o_rows));
  to_device(d, h, up, st);
  const Crops c = cl.at(d);
  label_chain(st, c, crows, cl.crop, words, cl.end - cl.nc, t);
  sig_lists(st, c, crows, cl.crop, t);
  t.stage("score", st, [&] {  // row partials and the terms
    for_each(st, rows, ScoreRows{ptr<Obj>(d, o_obj), c.pl, ptr<i64>(d, o_rows), n, words, ptr<i64>(d, o_parts)});
    for_each(st, n, ObjTerms{ptr<Obj>(d, o_obj), c, ptr<i64>(d, o_parts), ptr<ObjOut>(d, o_out)});
  });
  compact_sigs(st, c, ptr<i64>(d, o_ss), t);
  const size_t res = sizeof(ObjOut) * size_t(n), res_pad = (res + 255) / 256 * 256;
  u8* ho = L.out.get(res_pad + sizeof(i64) * size_t(n + 1));
  to_host(ho, d + o_out, res, st);
  to_host(ho + res_pad, d + o_ss, sizeof(i64) * size_t(n + 1), st);
  wait(st);
  const ObjOut* r = ptr<ObjOut>(ho, 0);
  for (i64 o = 0; o < n; ++o) {
    counters().label_pairs += r[o].relabelled;
    out.area.push_back(r[o].area), out.sum_x.push_back(r[o].sum_x), out.sum_y.push_back(r[o].sum_y);
    out.iou_tc.push_back(r[o].iou_tc), out.f5.push_back(r[o].f5), out.g_area.push_back(r[o].g_area);
    out.topo_freq.push_back(r[o].topo_freq);
  }
  out.sig_start.assign(ptr<i64>(ho + res_pad, 0), ptr<i64>(ho + res_pad, 0) + n + 1);
  out.sig = read_sigs(L, c, out.sig_start.back());
}

// ===================================================================== expert tables built on the device
// add_tables_from_raters (D7b): each cache's tables from its k rater planes, on the device, in one call for many
// caches: union, core and strict-majority vote word by word; the distance map from two exact EDTs (outside the union
// the distance to the union, inside the core the distance to the nearest non-core pixel, 0 in between; _BIG = 1e6
// when the EDT has no source, as objdel.reward.cache.edt); the perimeter max(1, 4-adjacent boundary pixels of the
// union); rater areas; the rater signatures from the pixel label chain on full-frame crops (so the crop layout needs
// no host round trip); voters (non-empty signature, or every rater when all are empty). Each cache's table image is
// built in the call's arena in its block layout (place) and copied into its own pool block after the call's one
// wait, when the signature lengths are known.

// One table set in one block: planes, map, per-rater arrays, the signatures last (their length is known last).
struct Placed {
  size_t size, r, u, c, v, d, a, ss, vt, s;
};

Placed place(i32 k, i32 H, i32 W, i64 nsig) {
  const i64 plane = i64(H) * words_per_row(W), px = i64(H) * W;
  Layout l;
  Placed q{};
  q.r = l.put<u32>(k * plane), q.u = l.put<u32>(plane), q.c = l.put<u32>(plane), q.v = l.put<u32>(plane);
  q.d = l.put<float>(px), q.a = l.put<i64>(k), q.ss = l.put<i64>(k + 1), q.vt = l.put<u8>(k), q.s = l.put<i32>(nsig);
  q.size = l.size;
  return q;
}

// Registers a filled block (called with the mutex held) -> its handle.
i64 register_tables(u8* b, const Placed& q, i32 H, i32 W, i32 k, i64 n_vote, i64 perimeter, double t5, i32 min_area) {
  Registry& r = registry();
  const TableView v{ptr<u32>(b, q.r), ptr<u32>(b, q.u), ptr<u32>(b, q.c),  ptr<u32>(b, q.v), ptr<float>(b, q.d),
                    ptr<i64>(b, q.a), ptr<i64>(b, q.ss), ptr<i32>(b, q.s), ptr<u8>(b, q.vt), n_vote, perimeter,
                    t5,               H,                W,                k,                min_area};
  const i64 handle = i64(r.block.size());
  r.block.push_back(b), r.bytes.push_back(q.size), r.view.push_back(v);
  ++r.live;
  counters().table_bytes += i64(q.size);
  return handle;
}

struct Build {  // one cache of the call: its uploaded rater planes, its table image and its scratch (device pointers)
  const u32* in;                       // k planes
  u32 *raters, *uni, *core, *vote;     // the image's planes
  float* dist;
  i64 *area, *ss;
  u8* voter;
  u32* notcore;                        // ~core inside the frame: the second EDT's sources
  const i32 *d_out, *d_in;             // the two squared EDTs (kNoSource: no source)
  i64* part;                           // H rows of k + 2 partials: rater areas, vote area, union boundary pixels
  i64 a0;                              // its first rater among the call's raters (= its first label-chain plane)
  i32 H, W, k;
};

HD inline u32 row_word(const u32* plane, i32 H, i32 wpr, i32 y, i32 j) {  // 0 outside the frame
  return y < 0 || y >= H || j < 0 || j >= wpr ? 0u : plane[i64(y) * wpr + j];
}

// Per call row: the image's rater, union, core and vote words, the non-core sources and the row partials. Vote bit:
// 2 x (raters with the bit) > k. Boundary: a union pixel with a 4-neighbour outside the union (outside the frame
// counts as outside), by word shifts.
struct TableRows {
  const Build* b;
  const i64* row_start;
  i64 n;
  HD u32 uni(const Build& t, i32 wpr, i32 y, i32 j) const {
    u32 u = 0;
    for (i32 r = 0; r < t.k; ++r) u |= row_word(t.in + r * i64(t.H) * wpr, t.H, wpr, y, j);
    return u;
  }
  HD void operator()(i64 i) const {
    const i64 o = segment(row_start, n, i);
    const Build& t = b[o];
    const i32 y = i32(i - row_start[o]), wpr = words_per_row(t.W);
    const i64 plane = i64(t.H) * wpr, at = i64(y) * wpr;
    i64* p = t.part + i64(y) * (t.k + 2);
    for (i32 r = 0; r < t.k + 2; ++r) p[r] = 0;
    u32 w[32];
    for (i32 j = 0; j < wpr; ++j) {
      u32 u = 0, c = ~0u, v = 0;
      for (i32 r = 0; r < t.k; ++r) {
        w[r] = t.in[r * plane + at + j];
        t.raters[r * plane + at + j] = w[r];
        u |= w[r], c &= w[r], p[r] += popc(w[r]);
      }
      for (u32 m = u; m; m &= m - 1) {
        const i32 bit = ctz(m);
        i32 cnt = 0;
        for (i32 r = 0; r < t.k; ++r) cnt += i32(w[r] >> bit & 1u);
        if (2 * cnt > t.k) v |= 1u << bit;
      }
      const i32 tail = t.W - 32 * j;
      const u32 valid = tail >= 32 ? ~0u : (1u << tail) - 1u;
      t.uni[at + j] = u, t.core[at + j] = c, t.vote[at + j] = v, t.notcore[at + j] = ~c & valid;
      const u32 left = (u << 1) | (uni(t, wpr, y, j - 1) >> 31), right = (u >> 1) | (uni(t, wpr, y, j + 1) << 31);
      p[t.k] += popc(v), p[t.k + 1] += popc(u & ~(uni(t, wpr, y - 1, j) & uni(t, wpr, y + 1, j) & left & right));
    }
  }
};

// The exact squared EDT of a frame's source bits (K2-draft EdtRowsH / EdtColsV): pass 1 along rows, pass 2 the lower
// envelope of parabolas down columns (Felzenszwalb and Huttenlocher 2012, Algorithm 1) in Meijster et al.'s integer
// form, floor division. Integers only, so any schedule gives the same values.
constexpr i32 kNoSource = 0x7FFFFFFF;  // squared distance when the frame has no source pixel

struct EdtFrame {
  const u32* src;    // H rows of words_per_row(W) words, bits past W are 0
  u16 *g, *s, *t;    // H x W scratch: row distances, parabola rows, their first rows (H + W <= 65534)
  i32* d;            // H x W squared distances
  i32 H, W;
};

// Pass 1, one warp per row: g = distance along the row to the nearest source, H + W when the row has none. Lane b
// takes pixel 32 j + b; the nearest source at or left of it is the highest set bit at or below b (clz), else the last
// source of the earlier words; the right side likewise (ctz). Each lane rewrites only its own pixels.
struct EdtRows {
  const EdtFrame* f;
  const i64* row_start;
  i64 n;
  HD void operator()(i64 i, Warp w) const {
    const i64 m = segment(row_start, n, i);
    const EdtFrame& e = f[m];
    const i32 y = i32(i - row_start[m]), W = e.W, wpr = words_per_row(W), inf = e.H + W;
    const u32* row = e.src + i64(y) * wpr;
    u16* o = e.g + i64(y) * W;
    i64 last = -1;  // column of the last source in the words before j
    for (i32 j = 0; j < wpr; ++j) {
      const u32 bits = row[j];
      for (i32 b = w.lane; b < 32; b += w.size) {
        const i32 x = 32 * j + b;
        if (x >= W) continue;
        const u32 low = bits & ((2u << b) - 1u);  // bits 0..b (2u << 31 wraps to 0: all bits)
        const i64 src = low ? i64(32 * j + 31 - clz(low)) : last;
        o[x] = u16(src >= 0 ? i32(x - src) : inf);
      }
      if (bits) last = 32 * j + 31 - clz(bits);
    }
    i64 next = -1;  // column of the first source in the words after j
    for (i32 j = wpr - 1; j >= 0; --j) {
      const u32 bits = row[j];
      for (i32 b = w.lane; b < 32; b += w.size) {
        const i32 x = 32 * j + b;
        if (x >= W) continue;
        const u32 high = bits & ~((1u << b) - 1u);  // bits b..31
        const i64 src = high ? i64(32 * j + ctz(high)) : next;
        if (src >= 0 && src - x < o[x]) o[x] = u16(src - x);
      }
      if (bits) next = 32 * j + ctz(bits);
    }
  }
};

// Pass 2, one thread per column: d(y) = min over rows i of (y - i)^2 + g(i)^2. Sep(a, u) = floor((u^2 - a^2 + g(u)^2 -
// g(a)^2) / (2 (u - a))) is the last row where parabola a is at least as low as parabola u (the numerator can be
// negative: floor, not truncation). A column of only H + W rows belongs to a frame with no source: kNoSource.
struct EdtCols {
  const EdtFrame* f;
  const i64* col_start;
  i64 n;
  HD void operator()(i64 c) const {
    const i64 m = segment(col_start, n, c);
    const EdtFrame& e = f[m];
    const i32 x = i32(c - col_start[m]), H = e.H, W = e.W;
    const i64 inf = H + W;
    auto G = [&](i64 r) -> i64 { return e.g[r * W + x]; };
    auto S = [&](i64 q) -> u16& { return e.s[q * W + x]; };
    auto T = [&](i64 q) -> u16& { return e.t[q * W + x]; };
    bool any = false;
    for (i32 r = 0; r < H && !any; ++r) any = G(r) < inf;
    if (!any) {
      for (i32 r = 0; r < H; ++r) e.d[i64(r) * W + x] = kNoSource;
      return;
    }
    auto F = [&](i64 r, i64 i) { const i64 gi = G(i); return (r - i) * (r - i) + gi * gi; };
    auto sep = [&](i64 a, i64 u) {
      const i64 ga = G(a), gu = G(u);
      return floordiv<i64>(u * u - a * a + gu * gu - ga * ga, 2 * (u - a));
    };
    i32 q = 0;
    S(0) = 0, T(0) = 0;
    for (i32 u = 1; u < H; ++u) {
      while (q >= 0 && F(T(q), S(q)) > F(T(q), u)) --q;
      if (q < 0) {
        q = 0, S(0) = u16(u);
      } else {
        const i64 v = 1 + sep(S(q), u);
        if (v < H) ++q, S(q) = u16(u), T(q) = u16(v);
      }
    }
    for (i32 u = H - 1; u >= 0; --u) {
      e.d[i64(u) * W + x] = i32(F(u, S(q)));
      if (u == T(q)) --q;
    }
  }
};

// float32(sqrt(float64(d2))): scipy's sqrt of the exact integer, rounded once more to float32 (cache.edt).
HD inline float edt_value(i32 d2) { return d2 == kNoSource ? 1e6f : float(sqrt(double(d2))); }

// Per pixel of one map: edt_value of its squared distance (gc.edt).
struct EdtMap {
  const i32* d2;
  float* out;
  HD void operator()(i64 p) const { out[p] = edt_value(d2[p]); }
};

// Per call row: the map where(union, dist_in_core, dist_out_union).
struct TableMap {
  const Build* b;
  const i64* row_start;
  i64 n;
  HD void operator()(i64 i) const {
    const i64 o = segment(row_start, n, i);
    const Build& t = b[o];
    const i32 y = i32(i - row_start[o]), wpr = words_per_row(t.W);
    for (i32 x = 0; x < t.W; ++x) {
      const i64 p = i64(y) * t.W + x, w = i64(y) * wpr + x / 32;
      const bool in_u = t.uni[w] >> (x % 32) & 1u, in_c = t.core[w] >> (x % 32) & 1u;
      t.dist[p] = in_u ? (in_c ? edt_value(t.d_in[p]) : 0.0f) : edt_value(t.d_out[p]);
    }
  }
};

// Per cache, after the label chain and compact_sigs (signatures one after the other, starts in ss): the rater areas,
// voters and signature starts into the image, and its scalars (rater areas, n_vote, perimeter, voters) into res.
struct TableTail {
  const Build* b;
  const i64* ss;
  i64 m, K;
  i64* res;    // K rater areas, then m vote areas, then m perimeters
  u8* voter;   // K
  HD void operator()(i64 o) const {
    const Build& t = b[o];
    i64 nv = 0, per = 0;
    bool any = false;
    for (i32 r = 0; r < t.k; ++r) {
      i64 a = 0;
      for (i32 y = 0; y < t.H; ++y) a += t.part[i64(y) * (t.k + 2) + r];
      t.area[r] = res[t.a0 + r] = a;
      t.ss[r] = ss[t.a0 + r] - ss[t.a0];
      any = any || ss[t.a0 + r + 1] > ss[t.a0 + r];
    }
    t.ss[t.k] = ss[t.a0 + t.k] - ss[t.a0];
    for (i32 y = 0; y < t.H; ++y) nv += t.part[i64(y) * (t.k + 2) + t.k], per += t.part[i64(y) * (t.k + 2) + t.k + 1];
    for (i32 r = 0; r < t.k; ++r) t.voter[r] = voter[t.a0 + r] = u8(!any || ss[t.a0 + r + 1] > ss[t.a0 + r]);
    res[K + o] = nv, res[K + m + o] = per > 1 ? per : 1;
  }
};

// ===================================================================== ladder certificate
// certify_polygons: candidate polygons rasterized as score_polygons does, compared with resident reference planes
// word by word, and every raster's size-filtered signature from the label chain; only counts and signatures leave.

// Per candidate frame row: an illegal shape's row is cleared (it draws nothing), then the popcounts of (c & r),
// (c | r) and (c ^ r) against its reference row are added to the shape's three counters (zeroed).
struct Overlap {
  Shapes s;
  const i64* defects;  // (n, kNumDefects)
  const i64* refw;     // per shape: the first word of its reference plane in `words`
  u32* cand;           // the candidates' rasters, at s.word_start
  const u32* words;
  i64* ov;
  HD void operator()(i64 i) const {
    i64 sh = 0;
    i32 row = 0;
    u32* m = cand + row_of(s, i, sh, row);
    const i32 wpr = words_per_row(s.W[sh]);
    bool bad = false;
    for (int k = 0; k < kNumDefects; ++k) bad = bad || defects[sh * kNumDefects + k] != 0;
    const u32* r = words + refw[sh] + i64(row) * wpr;
    i64 a = 0, b = 0, x = 0;
    for (i32 w = 0; w < wpr; ++w) {
      if (bad) m[w] = 0;
      a += popc(m[w] & r[w]), b += popc(m[w] | r[w]), x += popc(m[w] ^ r[w]);
    }
    atomic_add(ov + 3 * sh, a), atomic_add(ov + 3 * sh + 1, b), atomic_add(ov + 3 * sh + 2, x);
  }
};

struct Diff {  // one plane a & ~b of a pair query: first words of the result and of a and b, words per row
  i64 out, a, b;
  i32 wpr;
};

struct AndNot {  // per row of a pair plane
  const Diff* q;
  const i64* row_start;
  i64 n;
  u32* words;
  HD void operator()(i64 i) const {
    const i64 k = segment(row_start, n, i);
    const Diff& f = q[k];
    const i64 y = (i - row_start[k]) * f.wpr;
    for (i32 w = 0; w < f.wpr; ++w) words[f.out + y + w] = words[f.a + y + w] & ~words[f.b + y + w];
  }
};

// ===================================================================== topology
// topology: objdel.checker.topology on the label chain. After the filter and the final labels, the components and
// holes of every crop are ranked in raster order of their roots (first pixels) by a scan over crop rows.

// Per crop row: its foreground roots (components) and hole roots (bounded background components, not root 0).
struct RootCount {
  Crops c;
  i32 *comps, *holes;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Plane& b = c.pl[c.row(r, y)];
    const i32 cw = b.bw + 2;
    i32 nc = 0, nh = 0;
    for (i32 p = y * cw; p < (y + 1) * cw; ++p) {
      const i64 q = b.crop + p;
      nc += c.img[q] == 1 && c.fg[q] == p;
      nh += p != 0 && c.img[q] == 0 && c.bg[q] == p;
    }
    comps[r] = nc, holes[r] = nh;
  }
};

// Per crop row: each component's area (aux) and hole count (cnt) at its rank; cnt[root] becomes the rank.
struct CompEmit {
  Crops c;
  const i32* rank;  // exclusive scan of RootCount's comps
  i32 *cnt, *area, *nholes;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Plane& b = c.pl[c.row(r, y)];
    const i32 cw = b.bw + 2;
    i32 k = rank[r];
    for (i32 p = y * cw; p < (y + 1) * cw; ++p) {
      const i64 q = b.crop + p;
      if (c.img[q] != 1 || c.fg[q] != p) continue;
      area[k] = c.aux[q], nholes[k] = cnt[q], cnt[q] = k, ++k;
    }
  }
};

// Per crop row: each hole's area and the rank of the component above its root, at the hole's rank.
struct HoleEmit {
  Crops c;
  const i32* rank;  // exclusive scan of RootCount's holes
  const i32* cnt;   // component ranks at their roots (CompEmit)
  i32 *area, *comp;
  HD void operator()(i64 r) const {
    i32 y = 0;
    const Plane& b = c.pl[c.row(r, y)];
    const i32 cw = b.bw + 2;
    i32 h = rank[r];
    for (i32 p = y * cw; p < (y + 1) * cw; ++p) {
      const i64 q = b.crop + p;
      if (p == 0 || c.img[q] != 0 || c.bg[q] != p) continue;
      area[h] = c.aux[q], comp[h] = cnt[b.crop + c.fg[q - cw]], ++h;
    }
  }
};


// ===================================================================== evaluation sets
// eval_sets (D7): the pixel work of objdel.eval.metrics' item-level comparisons, for sets of K sample polygons against
// one cache's tables (k rater planes and the vote). Only integer counts come back (and the IoU_tc of each set's pixel
// vote, scored on the device by score_words): the M x M intersections of the drawn samples (raw raster AND legal),
// the raters and the vote; the K x K intersections of the RAW samples (duplicates); the M x M intersections of their
// boundary-IoU bands; per threshold the soft-Dice counts; per mask the surfels by length class on the (H + 1) x (W + 1)
// corner grid; per (sample, rater) the surfels farther than the NSD cut from the other surface, by an exact squared
// EDT of each surface (EdtRows, EdtCols). Every count is an integer sum, so any schedule and any chunk split give the
// same values. Per set, its plane table holds M drawn / rater / vote planes (an illegal sample's drawn plane is the
// chunk's zero plane), then K raw, M band and M eroded planes.

struct ESet {     // one set of a chunk (device pointers)
  i64 p0;         // its first entry in the plane table (3 M + K entries)
  u32* grid;      // (K + k) masks x 4 corner-grid planes: border, then the classes sqrt(2)/2, 1, sqrt(2)
  const i32* d2;  // (K + k) squared EDTs of the borders, (H + 1) x (W + 1) each
  u32* vote;      // the pixel-vote plane (lane 0's dev arena)
  i64 nsd_cut;
  i32 H, W, d, legal;  // frame, boundary-IoU band width, legal samples
};

HD inline u32 valid_bits(i32 W, i32 j) {  // the pixels of word j of a row of W pixels
  const i32 tail = W - 32 * j;
  return tail >= 32 ? ~0u : (1u << tail) - 1u;
}

// Per set row: the pixel vote (2 x legal drawn samples with the pixel > legal samples, objdel.eval.select.vote_mask)
// and, per threshold t, |P_t|, |R_t|, |P_t & R_t| with P_t = (drawn samples >= cut[t][0]), R_t = (raters >= cut[t][1]).
struct EvalPix {
  const ESet* s;
  u32* const* pp;
  const i64* start;
  i64 n;
  i32 K, k, T;
  const i32* cut;  // (n, T, 2)
  i64* soft;       // (n, T, 3)
  HD void operator()(i64 i) const {
    const i64 e = segment(start, n, i);
    const ESet& v = s[e];
    const i32 y = i32(i - start[e]), wpr = words_per_row(v.W);
    u32* const* p = pp + v.p0;
    const i32* c = cut + e * T * 2;
    i64 acc[3 * kMaxThresholds] = {};
    for (i32 j = 0; j < wpr; ++j) {
      const i64 at = i64(y) * wpr + j;
      const u32 valid = valid_bits(v.W, j);
      u32 vote = 0;
      for (i32 b = 0; b < 32 && (valid >> b & 1u); ++b) {
        i32 cs = 0, cr = 0;
        for (i32 q = 0; q < K; ++q) cs += i32(p[q][at] >> b & 1u);
        for (i32 q = 0; q < k; ++q) cr += i32(p[K + q][at] >> b & 1u);
        if (2 * cs > v.legal) vote |= 1u << b;
        for (i32 t = 0; t < T; ++t) {
          const bool pb = cs >= c[2 * t], rb = cr >= c[2 * t + 1];
          acc[3 * t] += pb, acc[3 * t + 1] += rb, acc[3 * t + 2] += pb && rb;
        }
      }
      v.vote[at] = vote;
    }
    for (i32 t = 0; t < 3 * T; ++t)
      if (acc[t]) atomic_add(soft + e * T * 3 + t, acc[t]);
  }
};

// Per (set, mask m of the M, row): the row erosion of m by d pixels (outside the frame counts as background), into its
// eroded plane. d <= 31, so a pixel's window lies in its word and the two next to it.
struct EvalErode {
  const ESet* s;
  u32* const* pp;
  const i64* start;
  i64 n;
  i32 K, M;
  HD void operator()(i64 i) const {
    const i64 e = segment(start, n, i);
    const ESet& v = s[e];
    const i64 local = i - start[e];
    const i32 m = i32(local / v.H), y = i32(local % v.H), wpr = words_per_row(v.W);
    const u32* a = pp[v.p0 + m] + i64(y) * wpr;
    u32* o = pp[v.p0 + 2 * M + K + m] + i64(y) * wpr;
    for (i32 j = 0; j < wpr; ++j) {
      const u32 w = a[j], l = j > 0 ? a[j - 1] : 0u, r = j + 1 < wpr ? a[j + 1] : 0u;
      u32 x = w;
      for (i32 dx = 1; dx <= v.d; ++dx) x &= ((w >> dx) | (r << (32 - dx))) & ((w << dx) | (l >> (32 - dx)));
      o[j] = x;
    }
  }
};

// Per (set, mask m, row): band = m & ~(the column erosion of m's eroded rows), metrics.band: the official
// mask_to_boundary (cv2.erode of the zero-padded mask, d iterations of a 3 x 3 square = one (2d + 1) square).
struct EvalBand {
  const ESet* s;
  u32* const* pp;
  const i64* start;
  i64 n;
  i32 K, M;
  HD void operator()(i64 i) const {
    const i64 e = segment(start, n, i);
    const ESet& v = s[e];
    const i64 local = i - start[e];
    const i32 m = i32(local / v.H), y = i32(local % v.H), wpr = words_per_row(v.W);
    const u32 *a = pp[v.p0 + m], *er = pp[v.p0 + 2 * M + K + m];
    u32* o = pp[v.p0 + M + K + m] + i64(y) * wpr;
    const bool inside = y - v.d >= 0 && y + v.d < v.H;
    for (i32 j = 0; j < wpr; ++j) {
      u32 x = inside ? ~0u : 0u;
      for (i32 dy = -v.d; inside && dy <= v.d; ++dy) x &= er[i64(y + dy) * wpr + j];
      o[j] = a[i64(y) * wpr + j] & ~x;
    }
  }
};

// One warp per (set, pair i <= j) of a group of G planes from entry `off` of each set's plane table:
// out[set][i][j] = out[set][j][i] = popcount(plane i & plane j) (K2-draft SetInter).
struct PairCount {
  const ESet* s;
  u32* const* pp;
  i64 pairs;  // per set: G (G + 1) / 2
  i32 off, G;
  i64* out;   // (n, G, G)
  HD void operator()(i64 t, Warp w) const {
    const i64 e = t / pairs;
    i64 local = t % pairs, i = 0;
    while (local >= G - i) local -= G - i, ++i;
    const i64 j = i + local;
    const ESet& v = s[e];
    const i64 nw = i64(v.H) * words_per_row(v.W);
    const u32 *a = pp[v.p0 + off + i], *b = pp[v.p0 + off + j];
    i64 c = 0;
    for (i64 q = w.lane; q < nw; q += w.size) c += popc(a[q] & b[q]);
    c = wsum(c);
    if (w.lane == 0) out[(e * G + i) * G + j] = c, out[(e * G + j) * G + i] = c;
  }
};

// Per (set, mask q of the K drawn samples and k raters, corner row y in [0, H]): the 2 x 2 code of every corner x in
// [0, W] (metrics._surface: 8 m(y-1, x-1) + 4 m(y-1, x) + 2 m(y, x-1) + m(y, x), outside the frame 0), bit-sliced: an odd
// count is a surfel of length sqrt(2)/2, two diagonal pixels (codes 6, 9) one of sqrt(2), two others one of 1; the
// border is every code but 0 and 15. Writes the four grid planes and adds the class counts to surf.
struct EvalGrid {
  const ESet* s;
  u32* const* pp;
  const i64* start;
  i64 n;
  i32 K, k;
  i64* surf;  // (n, K + k, 3)
  HD void operator()(i64 i) const {
    const i64 e = segment(start, n, i);
    const ESet& v = s[e];
    const i64 local = i - start[e];
    const i32 q = i32(local / (v.H + 1)), y = i32(local % (v.H + 1)), wpr = words_per_row(v.W),
              gwpr = words_per_row(v.W + 1);
    const i64 gp = i64(v.H + 1) * gwpr;
    const u32* m = pp[v.p0 + q];
    u32* g = v.grid + 4 * q * gp + i64(y) * gwpr;
    auto word = [&](i32 row, i32 j) { return row < 0 || row >= v.H || j < 0 || j >= wpr ? 0u : m[i64(row) * wpr + j]; };
    i64 n0 = 0, n1 = 0, n2 = 0;
    for (i32 j = 0; j < gwpr; ++j) {
      const u32 up = word(y - 1, j), dn = word(y, j);
      const u32 a = (up << 1) | (word(y - 1, j - 1) >> 31), b = up, c = (dn << 1) | (word(y, j - 1) >> 31), d = dn;
      const u32 odd = a ^ b ^ c ^ d, two = ~odd & (a | b | c | d) & ~(a & b & c & d), diag = two & ((a & d) | (b & c));
      g[0] = odd | two, g[gp] = odd, g[2 * gp] = two & ~diag, g[3 * gp] = diag;
      n0 += popc(odd), n1 += popc(two & ~diag), n2 += popc(diag);
      ++g;
    }
    i64* o = surf + (e * (K + k) + q) * 3;
    if (n0) atomic_add(o, n0);
    if (n1) atomic_add(o + 1, n1);
    if (n2) atomic_add(o + 2, n2);
  }
};

// Per (set, sample i, rater j, corner row y): the surfels of i farther than the cut from j's surface (direction 0) and
// those of j farther from i's (direction 1), by class. An empty surface (no source: kNoSource) is farther than any cut.
struct EvalNsd {
  const ESet* s;
  const i64* start;
  i64 n;
  i32 K, k;
  i64* nsd;  // (n, K, k, 2, 3)
  HD void operator()(i64 i) const {
    const i64 e = segment(start, n, i);
    const ESet& v = s[e];
    const i64 local = i - start[e];
    const i64 pair = local / (v.H + 1);
    const i32 a = i32(pair / k), b = K + i32(pair % k), y = i32(local % (v.H + 1)), gwpr = words_per_row(v.W + 1),
              gw = v.W + 1;
    const i64 gp = i64(v.H + 1) * gwpr, px = i64(v.H + 1) * gw;
    i64 acc[6] = {};
    for (int dir = 0; dir < 2; ++dir) {
      const i32 from = dir ? b : a, to = dir ? a : b;
      const i32* d2 = v.d2 + to * px + i64(y) * gw;
      for (int r = 0; r < 3; ++r) {
        const u32* g = v.grid + (4 * from + 1 + r) * gp + i64(y) * gwpr;
        for (i32 j = 0; j < gwpr; ++j)
          for (u32 m = g[j]; m; m &= m - 1) {
            const i32 d = d2[32 * j + ctz(m)];
            acc[3 * dir + r] += d == kNoSource || d > v.nsd_cut;
          }
      }
    }
    i64* o = nsd + ((e * K + a) * k + (b - K)) * 6;
    for (int r = 0; r < 6; ++r)
      if (acc[r]) atomic_add(o + r, acc[r]);
  }
};

// One chunk of eval_sets: sets [a, b) and their K b - K a shapes. The uploaded part first (the polygon inputs, the
// sets, the plane table, the cuts, the item starts, the EDT frames and their starts), then the work planes, the zero
// plane and the outputs (zeroed together).
struct EvalChunk {
  PolyChunk c;
  i64 nb, rows, erows, grows, nrows, frames, frows, fcols, maxpw, band_words, grid_words, px;
  size_t sets, pp, cut, rs, es, gs, ns, fr, frs, fcs, up, wide_r, wide_v, raw, band, er, grid, g, s, t, d2, zero,
      inter, raw_i, band_i, soft, surf, nsd, size;
};

EvalChunk eval_chunk(const Polygons& in, const std::vector<i64>& ws, const std::vector<TableView>& tv, i32 K, i32 k,
                     i32 T, i64 a, i64 b) {
  EvalChunk x{};
  const i32 M = K + k + 1, Q = K + k;
  x.c = plan(in, ws, a * K, b * K);
  x.nb = b - a;
  for (i64 e = a; e < b; ++e) {
    const TableView& v = tv[size_t(e)];
    const i64 pw = i64(v.H) * words_per_row(v.W), gp = i64(v.H + 1) * words_per_row(v.W + 1),
              px = i64(v.H + 1) * (v.W + 1);
    x.rows += v.H, x.erows += M * i64(v.H), x.grows += Q * i64(v.H + 1), x.nrows += i64(K) * k * (v.H + 1);
    x.frames += Q, x.frows += Q * i64(v.H + 1), x.fcols += Q * i64(v.W + 1);
    x.maxpw = mx(x.maxpw, pw), x.band_words += M * pw, x.grid_words += 4 * Q * gp, x.px += Q * px;
  }
  const i64 nb = x.nb, NP = 3 * M + K;
  Layout l;
  put_inputs(l, x.c);
  x.sets = l.put<ESet>(nb), x.pp = l.put<u32*>(nb * NP), x.cut = l.put<i32>(nb * T * 2);
  x.rs = l.put<i64>(nb + 1), x.es = l.put<i64>(nb + 1), x.gs = l.put<i64>(nb + 1), x.ns = l.put<i64>(nb + 1);
  x.fr = l.put<EdtFrame>(x.frames), x.frs = l.put<i64>(x.frames + 1), x.fcs = l.put<i64>(x.frames + 1);
  x.up = l.size;
  x.wide_r = l.put<u8>(x.c.n), x.wide_v = l.put<u8>(x.c.n), x.raw = l.put<u32>(x.c.nwords);
  x.band = l.put<u32>(x.band_words), x.er = l.put<u32>(x.band_words), x.grid = l.put<u32>(x.grid_words);
  x.g = l.put<u16>(x.px), x.s = l.put<u16>(x.px), x.t = l.put<u16>(x.px), x.d2 = l.put<i32>(x.px);
  x.zero = l.put<u32>(x.maxpw), x.inter = l.put<i64>(nb * M * M), x.raw_i = l.put<i64>(nb * K * K);
  x.band_i = l.put<i64>(nb * M * M), x.soft = l.put<i64>(nb * T * 3), x.surf = l.put<i64>(nb * Q * 3);
  x.nsd = l.put<i64>(nb * K * k * 6);
  x.size = l.size;
  return x;
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

// The pools' reserved bytes are read live, and only once the engine has a context (stats() must not open one: forked
// workers read it).
Stats stats() {
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Stats s = counters();
#ifdef __CUDACC__
  if (e.ready) {
    int dev = 0;
    cudaMemPool_t dp = nullptr;
    cuda_check(cudaGetDevice(&dev));
    cuda_check(cudaDeviceGetDefaultMemPool(&dp, dev));
    s.pool_reserved = reserved(pool()), s.default_pool_reserved = reserved(dp);
  }
#endif
  return s;
}

void reset_stats() {
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Stats& s = counters();
  const Stats kept = s;
  s = Stats{};
  s.arena_bytes = kept.arena_bytes, s.table_bytes = kept.table_bytes;
  s.live_device_bytes = kept.live_device_bytes;
}

PolygonResult check_polygons(const Polygons& in, Times* times) {
  GC_RANGE("check_polygons");
  const i64 n = in.n;
  PolygonResult out;
  out.word_start = word_starts(n, in.W, in.H);
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
  e.start("check_polygons", n);
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

MaskResult check_masks(const Masks& in, Times* times, const MaskOptions& opt) {
  GC_RANGE("check_masks");
  if (opt.seed_rows < 0 || opt.chunk_pixels < 1 || opt.chunk_pixels > kMaxChunkPixels || opt.max_scratch_bytes < 1)
    throw std::invalid_argument("check_masks needs seed_rows >= 0, 1 <= chunk_pixels <= 2^30, max_scratch_bytes >= 1");
  MaskResult out;
  out.ring_start = {0};
  out.pts_start = {0};
  if (in.n == 0) return out;
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  Timer t{times, {}};
  e.start("check_masks", in.n);
  // Chunks of at most chunk_pixels vertices, (W + 1)(H + 1) per mask (that also bounds rows and empty masks); a larger
  // mask goes alone. So every count inside a chunk fits an i32.
  for (i64 m0 = 0; m0 < in.n;) {
    i64 m1 = m0, cost = 0;
    while (m1 < in.n) {
      const i64 c = (i64(in.W[m1]) + 1) * (i64(in.H[m1]) + 1);
      if (m1 > m0 && cost + c > opt.chunk_pixels) break;
      cost += c, ++m1;
    }
    mask_range(in, m0, m1, opt, out, e.lanes[0], t);
    m0 = m1;
  }
  if (times) t.read();
  return out;
}

std::vector<std::string> test_hooks() {
  std::vector<std::string> h;
#ifdef GC_POISON_SCAN
  h.emplace_back("GC_POISON_SCAN");
#endif
#ifdef GC_ALIAS_SCAN_FOR_TEST
  h.emplace_back("GC_ALIAS_SCAN_FOR_TEST");
#endif
#ifdef K1_MUTATE_FG4
  h.emplace_back("K1_MUTATE_FG4");
#endif
  if (K1_ORDER != 0) h.push_back("K1_ORDER=" + std::to_string(K1_ORDER));
  return h;
}

// Each table set is one pool block (so drop_tables frees it alone), filled from one pageable host buffer (never the
// pinned staging: the number of cold caches changes per step); the batch makes one host wait.
std::vector<i64> add_tables_many(const std::vector<Tables>& ts) {
  GC_RANGE("add_tables");
  const size_t m = ts.size();
  if (m == 0) return {};
  std::vector<Placed> pl(m);
  std::vector<size_t> at(m + 1, 0);
  for (size_t k = 0; k < m; ++k) {
    pl[k] = place(ts[k].k, ts[k].H, ts[k].W, ts[k].sig_start[ts[k].k]);
    at[k + 1] = at[k] + (pl[k].size + 255) / 256 * 256;
  }
  std::vector<u8> host(at[m], 0);
  for (size_t k = 0; k < m; ++k) {
    const Tables& t = ts[k];
    const Placed& q = pl[k];
    u8* h = host.data() + at[k];
    const i64 plane = i64(t.H) * words_per_row(t.W), px = i64(t.H) * t.W;
    for (i32 r = 0; r < t.k; ++r) pack_plane(t.raters + r * px, t.H, t.W, ptr<u32>(h, q.r) + r * plane);
    pack_plane(t.uni, t.H, t.W, ptr<u32>(h, q.u));
    pack_plane(t.core, t.H, t.W, ptr<u32>(h, q.c));
    pack_plane(t.vote, t.H, t.W, ptr<u32>(h, q.v));
    std::copy_n(t.dist, px, ptr<float>(h, q.d));
    std::copy_n(t.rater_area, t.k, ptr<i64>(h, q.a));
    std::copy_n(t.sig_start, t.k + 1, ptr<i64>(h, q.ss));
    std::copy_n(t.voter, t.k, ptr<u8>(h, q.vt));
    std::copy_n(t.sig, t.sig_start[t.k], ptr<i32>(h, q.s));
  }
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  e.start("add_tables", i64(m));
  const Stream st = e.lanes[0].stream;
  std::vector<u8*> d(m, nullptr);
  try {
    for (size_t k = 0; k < m; ++k) {
      d[k] = static_cast<u8*>(alloc(true, pl[k].size, st));
      to_device(d[k], host.data() + at[k], pl[k].size, st);
    }
    wait(st);
  } catch (...) {
    for (size_t k = 0; k < m; ++k)
      if (d[k]) release(true, d[k], pl[k].size, st);
    throw;
  }
  std::vector<i64> handles(m);
  for (size_t k = 0; k < m; ++k) {
    const Tables& t = ts[k];
    handles[k] = register_tables(d[k], pl[k], t.H, t.W, t.k, t.n_vote, t.perimeter, t.t5, t.min_area);
  }
  return handles;
}

// One upload (the rater planes and the call's descriptors, from pageable memory), the kernels on lane 0 in its score
// arena, one wait; then the scalars to pageable memory (the device is idle) and each image into its own pool block in
// stream order. Pinned staging is never used (D5c: the cold-cache count changes per step).
BuiltTables add_tables_from_raters(const std::vector<Raters>& in) {
  GC_RANGE("add_tables_from_raters");
  const i64 m = i64(in.size());
  BuiltTables out;
  out.rater_area_start = out.sig_start = {0};
  if (m == 0) return out;
  i64 K = 0, rows = 0, cols = 0, words = 0, px = 0;
  for (const Raters& c : in) {
    const i64 plane = i64(c.H) * words_per_row(c.W);
    K += c.k, rows += c.H, cols += c.W, words += c.k * plane, px += i64(c.H) * c.W;
    out.rater_area_start.push_back(K);
  }
  // Layout: the uploaded part first (label-chain planes and crop rows, builds, row starts, EDT frames and their row
  // and column starts, rater words), then the images, scratch and results.
  std::vector<Plane> pl;
  std::vector<i64> crow_start{0};
  i64 crop = 0;
  for (const Raters& c : in)
    for (i32 r = 0; r < c.k; ++r) {
      pl.push_back(Plane{0, c.min_area, 0, crop, c.W, c.H, 0, 0, c.W, c.H});
      crop += i64(c.W + 2) * (c.H + 2);
      crow_start.push_back(crow_start.back() + c.H + 2);
    }
  CropLayout cl{K, crop};
  Layout l;
  cl.front(l);
  const size_t o_b = l.put<Build>(m), o_rs = l.put<i64>(m + 1), o_f = l.put<EdtFrame>(2 * m),
               o_frs = l.put<i64>(2 * m + 1), o_fcs = l.put<i64>(2 * m + 1), o_in = l.put<u32>(words), up = l.size;
  std::vector<Placed> img(size_t(m), Placed{});
  std::vector<size_t> o_img(static_cast<size_t>(m)), o_nc(o_img), o_part(o_img);
  for (i64 o = 0; o < m; ++o) {
    const Raters& c = in[size_t(o)];
    img[size_t(o)] = place(c.k, c.H, c.W, 0);
    o_img[size_t(o)] = l.put<u8>(i64(img[size_t(o)].s));
  }
  const size_t img_end = l.size;
  for (i64 o = 0; o < m; ++o) {
    const Raters& c = in[size_t(o)];
    o_nc[size_t(o)] = l.put<u32>(i64(c.H) * words_per_row(c.W));
    o_part[size_t(o)] = l.put<i64>(i64(c.H) * (c.k + 2));
  }
  const size_t o_g = l.put<u16>(2 * px), o_s = l.put<u16>(2 * px), o_t = l.put<u16>(2 * px), o_d = l.put<i32>(2 * px);
  cl.back(l);
  const size_t o_ss = l.put<i64>(K + 1), o_res = l.put<i64>(K + 2 * m), o_vt = l.put<u8>(K);

  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  Timer t{nullptr, {}};
  e.start("add_tables_from_raters", m);
  Lane& L = e.lanes[0];
  const Stream st = L.stream;
  u8* d = L.score.get(l.size);
  std::vector<u8> host(up, 0);
  Build* bs = ptr<Build>(host.data(), o_b);
  i64* rs = ptr<i64>(host.data(), o_rs);
  EdtFrame* fs = ptr<EdtFrame>(host.data(), o_f);
  i64 *frs = ptr<i64>(host.data(), o_frs), *fcs = ptr<i64>(host.data(), o_fcs);
  rs[0] = frs[0] = fcs[0] = 0;
  i64 w = 0, pix = 0, a0 = 0;
  for (i64 o = 0; o < m; ++o) {
    const Raters& c = in[size_t(o)];
    const Placed& q = img[size_t(o)];
    const i64 plane = i64(c.H) * words_per_row(c.W), n = i64(c.H) * c.W;
    for (i32 r = 0; r < c.k; ++r) {
      pack_plane(c.px + r * n, c.H, c.W, ptr<u32>(host.data(), o_in) + w + r * plane);
      pl[size_t(a0 + r)].words = i64(o_in / sizeof(u32)) + w + r * plane;
    }
    u8* im = d + o_img[size_t(o)];
    const i32* d_out = ptr<i32>(d, o_d) + 2 * pix;
    bs[o] = Build{ptr<u32>(d, o_in) + w, ptr<u32>(im, q.r), ptr<u32>(im, q.u), ptr<u32>(im, q.c), ptr<u32>(im, q.v),
                  ptr<float>(im, q.d), ptr<i64>(im, q.a), ptr<i64>(im, q.ss), ptr<u8>(im, q.vt), ptr<u32>(d, o_nc[size_t(o)]),
                  d_out, d_out + n, ptr<i64>(d, o_part[size_t(o)]), a0, c.H, c.W, c.k};
    for (int f = 0; f < 2; ++f) {  // frame 2o: sources = the union (distance to it); 2o + 1: the non-core pixels
      const i64 at = 2 * pix + f * n;
      fs[2 * o + f] = EdtFrame{f ? bs[o].notcore : bs[o].uni, ptr<u16>(d, o_g) + at, ptr<u16>(d, o_s) + at,
                               ptr<u16>(d, o_t) + at, ptr<i32>(d, o_d) + at, c.H, c.W};
      frs[2 * o + f + 1] = frs[2 * o + f] + c.H, fcs[2 * o + f + 1] = fcs[2 * o + f] + c.W;
    }
    rs[o + 1] = rs[o] + c.H;
    w += c.k * plane, pix += n, a0 += c.k;
  }
  cl.stage(host.data(), pl, crow_start);
  to_device(d, host.data(), up, st);
  zero(d + o_img[0], img_end - o_img[0], st);  // the images' padding is copied into the blocks
  const Build* db = ptr<Build>(d, o_b);
  const i64* drs = ptr<i64>(d, o_rs);
  const EdtFrame* df = ptr<EdtFrame>(d, o_f);
  t.stage("tables", st, [&] {
    for_each(st, rows, TableRows{db, drs, m});
    for_each_warp(st, 2 * rows, EdtRows{df, ptr<i64>(d, o_frs), 2 * m});
    for_each(st, 2 * cols, EdtCols{df, ptr<i64>(d, o_fcs), 2 * m});
    for_each(st, rows, TableMap{db, drs, m});
  });
  const Crops c = cl.at(d);
  label_chain(st, c, crow_start.back(), crop, ptr<u32>(d, 0), cl.end - cl.nc, t);
  sig_lists(st, c, crow_start.back(), crop, t);
  compact_sigs(st, c, ptr<i64>(d, o_ss), t);
  t.stage("tables", st, [&] { for_each(st, m, TableTail{db, ptr<i64>(d, o_ss), m, K, ptr<i64>(d, o_res), ptr<u8>(d, o_vt)}); });
  wait(st);
  std::vector<i64> res(size_t(K + 2 * m));
  out.sig_start.resize(size_t(K + 1));
  out.voter.resize(size_t(K));
  read_back(res.data(), d + o_res, sizeof(i64) * res.size(), st);
  read_back(out.sig_start.data(), d + o_ss, sizeof(i64) * size_t(K + 1), st);
  read_back(out.voter.data(), d + o_vt, size_t(K), st);
  out.sig.resize(size_t(out.sig_start.back()));
  read_back(out.sig.data(), c.aux, sizeof(i32) * out.sig.size(), st);
  out.rater_area.assign(res.begin(), res.begin() + K);
  out.n_vote.assign(res.begin() + K, res.begin() + K + m);
  out.perimeter.assign(res.begin() + K + m, res.end());
  std::vector<u8*> blk(size_t(m), nullptr);
  std::vector<Placed> q(img);
  try {
    for (i64 o = 0; o < m; ++o) {
      const Raters& r = in[size_t(o)];
      const i64 s0 = out.sig_start[size_t(out.rater_area_start[size_t(o)])];
      const i64 ns = out.sig_start[size_t(out.rater_area_start[size_t(o + 1)])] - s0;
      q[size_t(o)] = place(r.k, r.H, r.W, ns);
      blk[size_t(o)] = static_cast<u8*>(alloc(true, q[size_t(o)].size, st));
      copy_device(blk[size_t(o)], d + o_img[size_t(o)], q[size_t(o)].s, st);
      copy_device(blk[size_t(o)] + q[size_t(o)].s, c.aux + s0, sizeof(i32) * size_t(ns), st);
    }
  } catch (...) {
    for (i64 o = 0; o < m; ++o)
      if (blk[size_t(o)]) release(true, blk[size_t(o)], q[size_t(o)].size, st);
    throw;
  }
  for (i64 o = 0; o < m; ++o) {
    const Raters& r = in[size_t(o)];
    out.handles.push_back(register_tables(blk[size_t(o)], q[size_t(o)], r.H, r.W, r.k, out.n_vote[size_t(o)],
                                          out.perimeter[size_t(o)], r.t5, r.min_area));
  }
  return out;
}

i64 add_tables(const Tables& t) { return add_tables_many({t})[0]; }

// One map's EDT with the table build's kernels (EdtRows, EdtCols, edt_value): the sources are the zero pixels. One
// upload from pageable memory, one wait, the map read back (the device is idle), as add_tables_from_raters.
std::vector<float> edt(const u8* nonzero, i32 H, i32 W) {
  GC_RANGE("edt");
  const i64 n = i64(H) * W, wpr = words_per_row(W);
  std::vector<float> out(static_cast<size_t>(n));
  if (n == 0) return out;
  Layout l;
  const size_t o_f = l.put<EdtFrame>(1), o_rs = l.put<i64>(2), o_cs = l.put<i64>(2), o_src = l.put<u32>(H * wpr),
               up = l.size;
  const size_t o_g = l.put<u16>(n), o_s = l.put<u16>(n), o_t = l.put<u16>(n), o_d = l.put<i32>(n),
               o_out = l.put<float>(n);
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  e.start("edt", 1);
  Lane& L = e.lanes[0];
  const Stream st = L.stream;
  u8* d = L.score.get(l.size);
  std::vector<u8> host(up, 0);
  u32* src = ptr<u32>(host.data(), o_src);
  for (i32 y = 0; y < H; ++y)
    for (i32 x = 0; x < W; ++x)
      if (!nonzero[i64(y) * W + x]) src[y * wpr + x / 32] |= u32(1) << (x % 32);
  *ptr<EdtFrame>(host.data(), o_f) = EdtFrame{ptr<u32>(d, o_src), ptr<u16>(d, o_g), ptr<u16>(d, o_s),
                                              ptr<u16>(d, o_t), ptr<i32>(d, o_d), H, W};
  ptr<i64>(host.data(), o_rs)[0] = ptr<i64>(host.data(), o_cs)[0] = 0;
  ptr<i64>(host.data(), o_rs)[1] = H, ptr<i64>(host.data(), o_cs)[1] = W;
  to_device(d, host.data(), up, st);
  const EdtFrame* df = ptr<EdtFrame>(d, o_f);
  for_each_warp(st, H, EdtRows{df, ptr<i64>(d, o_rs), 1});
  for_each(st, W, EdtCols{df, ptr<i64>(d, o_cs), 1});
  for_each(st, n, EdtMap{ptr<i32>(d, o_d), ptr<float>(d, o_out)});
  wait(st);
  read_back(out.data(), d + o_out, sizeof(float) * size_t(n), st);
  return out;
}

// Freed in stream order on lane 0 (no host wait); every call drains both lanes before it returns, so no kernel still
// reads the tables.
void drop_tables(i64 handle) {
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Registry& r = registry();
  if (handle < 0 || handle >= i64(r.block.size()) || !r.block[size_t(handle)])
    throw std::invalid_argument("unknown tables handle");
  clear_error();
  release(true, r.block[size_t(handle)], r.bytes[size_t(handle)], e.lanes[0].stream);
  counters().table_bytes -= i64(r.bytes[size_t(handle)]);
  r.block[size_t(handle)] = nullptr;
  --r.live;
}

void trim(i64 bytes) {
  if (bytes < 0) throw std::invalid_argument("trim needs bytes >= 0");
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
#ifdef __CUDACC__
  if (e.ready) cuda_check(cudaMemPoolTrimTo(pool(), size_t(bytes)));
#endif
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
  const std::vector<i64> word_start = word_starts(n, in.W, in.H);
  if (n == 0) return out;
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Drain drain{e};
  const std::vector<TableView> tv = views(tables, n, in.W, in.H);
  clear_error();
  Timer t{times, {}};
  e.start("score_polygons", n);
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
  const std::vector<i64> word_start = word_starts(n, in.W, in.H);
  if (n == 0) return out;
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Drain drain{e};
  const std::vector<TableView> tv = views(tables, n, in.W, in.H);
  clear_error();
  Timer t{times, {}};
  e.start("score_masks", n);
  Lane& L = e.lanes[0];
  const u32* words = upload_masks(L, in, word_start);
  out.defects.assign(size_t(n * kNumDefects), 0);
  score_words(L, tv, words, word_start, out, t);
  if (times) t.read();
  return out;
}

Certified certify_polygons(const Polygons& in, const Masks& ref, const i32* ref_of, const i64* min_area,
                           const Pairs& pr) {
  GC_RANGE("certify_polygons");
  const i64 n = in.n, P = pr.n;
  Certified out;
  out.sig_start = out.ref_sig_start = out.pair_sig_start = out.pair_ref_sig_start = {0};
  if (n == 0) return out;
  const std::vector<i64> word_start = word_starts(n, in.W, in.H), rws = word_starts(ref.n, ref.W, ref.H);
  std::vector<i64> pws(size_t(2 * P + 1), 0);
  std::vector<i32> qW(size_t(2 * P)), qH(size_t(2 * P));  // pair planes: P candidate pairs, then P reference pairs
  for (i64 q = 0; q < 2 * P; ++q) {
    const i64 a = q < P ? pr.shapes[2 * q] : pr.refs[2 * (q - P)];
    qW[size_t(q)] = q < P ? in.W[a] : ref.W[a], qH[size_t(q)] = q < P ? in.H[a] : ref.H[a];
    pws[size_t(q + 1)] = pws[size_t(q)] + i64(qH[size_t(q)]) * words_per_row(qW[size_t(q)]);
  }
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  Timer t{nullptr, {}};
  e.start("certify_polygons", n);
  Lane& L = e.lanes[0];
  const Stream st = L.stream;
  const PolyChunk c = plan(in, word_start, 0, n);
  // After the chunk in the checker arena: the reference planes, per shape its reference's first word, the pair
  // queries (copied in), the pair planes and the overlap counters. Words are indexed from the arena's start.
  Layout x;
  const size_t x_rw = x.put<u32>(rws.back()), x_refw = x.put<i64>(n), x_diff = x.put<Diff>(2 * P),
               x_drows = x.put<i64>(2 * P + 1), x_up = x.size, x_pw = x.put<u32>(pws.back()), x_ov = x.put<i64>(3 * n);
  const size_t base = (c.size + 255) / 256 * 256, xin = (c.in_bytes + 255) / 256 * 256;
  u8* d = L.dev.get(base + x.size);
  u8* h = L.in.get(xin + x_up);
  L.out.get(c.words - c.defects);
  pack(in, word_start, c, h);
  u8* hx = h + xin;
  const i64 rw0 = i64(base + x_rw) / 4, pw0 = i64(base + x_pw) / 4, cw0 = i64(c.words) / 4;
  std::fill_n(hx + x_rw, sizeof(u32) * size_t(rws.back()), u8(0));
  for (i64 j = 0; j < ref.n; ++j)
    pack_plane(ref.bits + ref.start[j], ref.H[j], ref.W[j], ptr<u32>(hx, x_rw) + rws[size_t(j)]);
  for (i64 s = 0; s < n; ++s) ptr<i64>(hx, x_refw)[s] = rw0 + rws[size_t(ref_of[s])];
  i64* drows = ptr<i64>(hx, x_drows);
  drows[0] = 0;
  for (i64 q = 0; q < 2 * P; ++q) {
    const i64 p = q % P;
    const i64 a = q < P ? cw0 + word_start[size_t(pr.shapes[2 * p])] : rw0 + rws[size_t(pr.refs[2 * p])];
    const i64 b = q < P ? cw0 + word_start[size_t(pr.shapes[2 * p + 1])] : rw0 + rws[size_t(pr.refs[2 * p + 1])];
    ptr<Diff>(hx, x_diff)[q] = Diff{pw0 + pws[size_t(q)], a, b, words_per_row(qW[size_t(q)])};
    drows[q + 1] = drows[q] + qH[size_t(q)];
  }
  launch(c, in.G, L, false, t);
  to_device(d + base, hx, x_up, st);
  zero(d + base + x_ov, sizeof(i64) * size_t(3 * n), st);
  u32* words = ptr<u32>(d, 0);
  t.stage("certify", st, [&] {
    for_each(st, c.rows, Overlap{shapes(c, d, in.G), ptr<i64>(d, c.defects), ptr<i64>(d, base + x_refw),
                                 ptr<u32>(d, c.words), words, ptr<i64>(d, base + x_ov)});
    for_each(st, drows[2 * P], AndNot{ptr<Diff>(d, base + x_diff), ptr<i64>(d, base + x_drows), 2 * P, words});
  });
  wait(st);
  for (i64 s = 0; s < n; ++s) {
    const i64* df = ptr<i64>(L.out.p, 0) + s * kNumDefects;
    out.legal.push_back(std::all_of(df, df + kNumDefects, [](i64 v) { return v == 0; }));
  }

  // The planes of the label chain: the candidates, each (reference, min_area) once, the pair planes.
  std::vector<Plane> pl;
  for (i64 s = 0; s < n; ++s) pl.push_back(Plane{cw0 + word_start[size_t(s)], min_area[s], 0, 0, in.W[s], in.H[s]});
  std::vector<i64> ref_plane(static_cast<size_t>(n));
  for (i64 s = 0; s < n; ++s) {
    i64 k = n;
    while (k < i64(pl.size()) && !(pl[size_t(k)].words == rw0 + rws[size_t(ref_of[s])] &&
                                   pl[size_t(k)].min_area == min_area[s]))
      ++k;
    if (k == i64(pl.size()))
      pl.push_back(Plane{rw0 + rws[size_t(ref_of[s])], min_area[s], 0, 0, ref.W[ref_of[s]], ref.H[ref_of[s]]});
    ref_plane[size_t(s)] = k;
  }
  const i64 q0 = i64(pl.size());
  for (i64 q = 0; q < 2 * P; ++q)
    pl.push_back(Plane{pw0 + pws[size_t(q)], pr.min_area[q % P], 0, 0, qW[size_t(q)], qH[size_t(q)]});
  const i64 np = i64(pl.size());
  std::vector<i64> crow_start;
  CropLayout cl = boxes(L, pl, words, crow_start, t);
  Layout l;
  cl.front(l);
  const size_t up = l.size;
  cl.back(l);
  const size_t o_ss = l.put<i64>(np + 1);
  u8* ds = L.score.get(l.size);
  cl.stage(L.in.get(up), pl, crow_start);
  to_device(ds, L.in.p, up, st);
  const Crops cr = cl.at(ds);
  label_chain(st, cr, crow_start.back(), cl.crop, words, cl.end - cl.nc, t);
  sig_lists(st, cr, crow_start.back(), cl.crop, t);
  compact_sigs(st, cr, ptr<i64>(ds, o_ss), t);
  const size_t ov = sizeof(i64) * size_t(3 * n), ov_pad = (ov + 255) / 256 * 256;
  u8* ho = L.out.get(ov_pad + sizeof(i64) * size_t(np + 1));
  to_host(ho, d + base + x_ov, ov, st);
  to_host(ho + ov_pad, ds + o_ss, sizeof(i64) * size_t(np + 1), st);
  wait(st);
  const i64* cnt = ptr<i64>(ho, 0);
  for (i64 s = 0; s < n; ++s)
    out.inter.push_back(cnt[3 * s]), out.uni.push_back(cnt[3 * s + 1]), out.equal.push_back(cnt[3 * s + 2] == 0);
  const std::vector<i64> ss(ptr<i64>(ho + ov_pad, 0), ptr<i64>(ho + ov_pad, 0) + np + 1);
  const std::vector<i32> sig = read_sigs(L, cr, ss.back());
  auto put = [&](i64 k, std::vector<i64>& start, std::vector<i32>& v) {
    v.insert(v.end(), sig.begin() + ss[size_t(k)], sig.begin() + ss[size_t(k + 1)]);
    start.push_back(i64(v.size()));
  };
  for (i64 s = 0; s < n; ++s) put(s, out.sig_start, out.sig), put(ref_plane[size_t(s)], out.ref_sig_start, out.ref_sig);
  for (i64 p = 0; p < P; ++p)
    put(q0 + p, out.pair_sig_start, out.pair_sig), put(q0 + P + p, out.pair_ref_sig_start, out.pair_ref_sig);
  return out;
}


// Sets in chunks whose scratch (one lane-0 score arena block) stays within the cap; per chunk one upload from pageable
// memory, the kernels, one wait and the counts read back (the device is idle). Then one score_words call scores every
// set's pixel vote (in lane 0's dev arena) against its tables.
EvalCounts eval_sets(const Polygons& in, const EvalSets& es) {
  GC_RANGE("eval_sets");
  const i64 n = es.n;
  const i32 K = es.K, T = es.T;
  EvalCounts out;
  if (n == 0) return out;
  const std::vector<i64> ws = word_starts(in.n, in.W, in.H);
  std::vector<i64> shape_tables(static_cast<size_t>(in.n)), vws(size_t(n + 1), 0);
  for (i64 q = 0; q < in.n; ++q) shape_tables[size_t(q)] = es.tables[q / K];
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  Drain drain{e};
  const std::vector<TableView> all = views(shape_tables.data(), in.n, in.W, in.H);
  const i32 k = out.k = all[0].k, M = K + k + 1, Q = K + k;
  std::vector<TableView> tv;
  for (i64 s = 0; s < n; ++s) {
    tv.push_back(all[size_t(s * K)]);
    if (tv.back().k != k) throw std::invalid_argument("eval_sets: every set's tables must have the same raters");
    for (i32 c = 0; c < T; ++c)
      if (es.cuts[(s * T + c) * 2 + 1] > k) throw std::invalid_argument("eval_sets: rater cuts must be in [0, k]");
    vws[size_t(s + 1)] = vws[size_t(s)] + i64(tv.back().H) * words_per_row(tv.back().W);
  }
  out.inter.assign(size_t(n * M * M), 0), out.inter_raw.assign(size_t(n * K * K), 0);
  out.band_inter.assign(size_t(n * M * M), 0), out.soft.assign(size_t(n * T * 3), 0);
  out.surf.assign(size_t(n * Q * 3), 0), out.nsd.assign(size_t(n * K * k * 6), 0);
  clear_error();
  Timer t{nullptr, {}};
  e.start("eval_sets", n);
  Lane& L = e.lanes[0];
  const Stream st = L.stream;
  u32* vote = ptr<u32>(L.dev.get(sizeof(u32) * size_t(vws[size_t(n)])), 0);
  const i64 NP = 3 * M + K;
  for (i64 a = 0, b = 0; a < n; a = b) {
    b = a + 1;
    if (i64(eval_chunk(in, ws, tv, K, k, T, a, b).size) > es.max_scratch_bytes)
      throw std::invalid_argument("eval_sets: one set needs more scratch than max_scratch_bytes");
    while (b < n && i64(eval_chunk(in, ws, tv, K, k, T, a, b + 1).size) <= es.max_scratch_bytes) ++b;
    const EvalChunk x = eval_chunk(in, ws, tv, K, k, T, a, b);
    const i64 nb = x.nb;
    counters().eval_scratch_peak = mx(counters().eval_scratch_peak, i64(x.size));
    u8* d = L.score.get(x.size);
    std::vector<u8> host(x.up, 0);
    pack(in, ws, x.c, host.data());
    ESet* sets = ptr<ESet>(host.data(), x.sets);
    u32** pp = ptr<u32*>(host.data(), x.pp);
    i64 *rs = ptr<i64>(host.data(), x.rs), *ers = ptr<i64>(host.data(), x.es), *grs = ptr<i64>(host.data(), x.gs),
        *nrs = ptr<i64>(host.data(), x.ns), *frs = ptr<i64>(host.data(), x.frs), *fcs = ptr<i64>(host.data(), x.fcs);
    EdtFrame* fr = ptr<EdtFrame>(host.data(), x.fr);
    rs[0] = ers[0] = grs[0] = nrs[0] = frs[0] = fcs[0] = 0;
    i64 band = 0, grid = 0, px = 0, f = 0;
    for (i64 j = 0; j < nb; ++j) {
      const i64 s = a + j;
      const TableView& v = tv[size_t(s)];
      const i64 pw = i64(v.H) * words_per_row(v.W), gp = i64(v.H + 1) * words_per_row(v.W + 1),
                fp = i64(v.H + 1) * (v.W + 1);
      u32** p = pp + j * NP;
      i32 legal = 0;
      for (i32 q = 0; q < K; ++q) {
        u32* raw = ptr<u32>(d, x.raw) + ws[size_t(s * K + q)] - ws[size_t(a * K)];
        const bool ok = es.legal[s * K + q] != 0;
        legal += ok;
        p[q] = ok ? raw : ptr<u32>(d, x.zero), p[M + q] = raw;
      }
      for (i32 r = 0; r < k; ++r) p[K + r] = const_cast<u32*>(v.raters) + r * pw;
      p[K + k] = const_cast<u32*>(v.vote);
      for (i32 m = 0; m < M; ++m)
        p[M + K + m] = ptr<u32>(d, x.band) + band + m * pw, p[2 * M + K + m] = ptr<u32>(d, x.er) + band + m * pw;
      for (i32 q = 0; q < Q; ++q) {
        fr[f] = EdtFrame{ptr<u32>(d, x.grid) + grid + 4 * q * gp, ptr<u16>(d, x.g) + px + q * fp,
                         ptr<u16>(d, x.s) + px + q * fp,          ptr<u16>(d, x.t) + px + q * fp,
                         ptr<i32>(d, x.d2) + px + q * fp,         v.H + 1,
                         v.W + 1};
        frs[f + 1] = frs[f] + v.H + 1, fcs[f + 1] = fcs[f] + v.W + 1, ++f;
      }
      sets[j] = ESet{j * NP, ptr<u32>(d, x.grid) + grid, ptr<i32>(d, x.d2) + px, vote + vws[size_t(s)],
                     es.nsd_cut[s], v.H, v.W, es.band_d[s], legal};
      std::copy_n(es.cuts + s * T * 2, T * 2, ptr<i32>(host.data(), x.cut) + j * T * 2);
      rs[j + 1] = rs[j] + v.H, ers[j + 1] = ers[j] + M * i64(v.H), grs[j + 1] = grs[j] + Q * i64(v.H + 1);
      nrs[j + 1] = nrs[j] + i64(K) * k * (v.H + 1);
      band += M * pw, grid += 4 * Q * gp, px += Q * fp;
    }
    to_device(d, host.data(), x.up, st);
    zero(d + x.zero, x.size - x.zero, st);
    const Shapes sh = shapes(x.c, d, in.G);
    const ESet* ds = ptr<ESet>(d, x.sets);
    u32* const* dp = ptr<u32*>(d, x.pp);
    const EdtFrame* df = ptr<EdtFrame>(d, x.fr);
    i64 *inter = ptr<i64>(d, x.inter), *raw_i = ptr<i64>(d, x.raw_i), *band_i = ptr<i64>(d, x.band_i);
    for_each(st, x.c.n, Precision{sh, ptr<u8>(d, x.wide_r), ptr<u8>(d, x.wide_v)});
    raster(st, sh, ptr<u8>(d, x.wide_r), ptr<u32>(d, x.raw), x.c.rows);
    for_each(st, x.rows, EvalPix{ds, dp, ptr<i64>(d, x.rs), nb, K, k, T, ptr<i32>(d, x.cut), ptr<i64>(d, x.soft)});
    for_each(st, x.erows, EvalErode{ds, dp, ptr<i64>(d, x.es), nb, K, M});
    for_each(st, x.erows, EvalBand{ds, dp, ptr<i64>(d, x.es), nb, K, M});
    const i64 pm = i64(M) * (M + 1) / 2, pk = i64(K) * (K + 1) / 2;
    for_each_warp(st, nb * pm, PairCount{ds, dp, pm, 0, M, inter});
    for_each_warp(st, nb * pk, PairCount{ds, dp, pk, M, K, raw_i});
    for_each_warp(st, nb * pm, PairCount{ds, dp, pm, M + K, M, band_i});
    for_each(st, x.grows, EvalGrid{ds, dp, ptr<i64>(d, x.gs), nb, K, k, ptr<i64>(d, x.surf)});
    for_each_warp(st, x.frows, EdtRows{df, ptr<i64>(d, x.frs), x.frames});
    for_each(st, x.fcols, EdtCols{df, ptr<i64>(d, x.fcs), x.frames});
    for_each(st, x.nrows, EvalNsd{ds, ptr<i64>(d, x.ns), nb, K, k, ptr<i64>(d, x.nsd)});
    wait(st);
    read_back(out.inter.data() + a * M * M, inter, sizeof(i64) * size_t(nb * M * M), st);
    read_back(out.inter_raw.data() + a * K * K, raw_i, sizeof(i64) * size_t(nb * K * K), st);
    read_back(out.band_inter.data() + a * M * M, band_i, sizeof(i64) * size_t(nb * M * M), st);
    read_back(out.soft.data() + a * T * 3, d + x.soft, sizeof(i64) * size_t(nb * T * 3), st);
    read_back(out.surf.data() + a * Q * 3, d + x.surf, sizeof(i64) * size_t(nb * Q * 3), st);
    read_back(out.nsd.data() + a * K * k * 6, d + x.nsd, sizeof(i64) * size_t(nb * K * k * 6), st);
  }
  Scores sc;
  score_words(L, tv, vote, vws, sc, t);
  out.vote_iou_tc = std::move(sc.iou_tc);
  return out;
}

Topology topology(const Masks& in, const i64* min_area, bool want_masks) {
  GC_RANGE("topology");
  const i64 n = in.n;
  Topology out;
  out.comp_start = out.hole_start = out.sig_start = {0};
  if (want_masks) out.masks.assign(size_t(in.start[n]), 0);
  if (n == 0) return out;
  const std::vector<i64> word_start = word_starts(n, in.W, in.H);
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  Timer t{nullptr, {}};
  e.start("topology", in.n);
  Lane& L = e.lanes[0];
  const Stream st = L.stream;
  const u32* words = upload_masks(L, in, word_start);
  std::vector<Plane> pl;
  for (i64 k = 0; k < n; ++k) pl.push_back(Plane{word_start[size_t(k)], min_area[k], 0, 0, in.W[k], in.H[k]});
  std::vector<i64> crow_start;
  CropLayout cl = boxes(L, pl, words, crow_start, t);
  if (cl.crop > INT_MAX) throw std::invalid_argument("topology: more than 2^31 - 1 crop pixels in one call");
  const i64 crows = crow_start.back();
  i64 fgpx = 0;
  for (const Plane& b : pl) fgpx += b.area;  // bounds the components; the other crop pixels bound the holes
  Layout l;
  cl.front(l);
  const size_t up = l.size;
  cl.back(l);
  const size_t o_cnt = l.put<i32>(cl.crop), o_rc = l.put<i32>(crows + 1), o_rh = l.put<i32>(crows + 1),
               o_kc = l.put<i32>(crows + 1), o_kh = l.put<i32>(crows + 1), o_ca = l.put<i32>(fgpx),
               o_ch = l.put<i32>(fgpx), o_ha = l.put<i32>(cl.crop - fgpx), o_hc = l.put<i32>(cl.crop - fgpx);
  u8* d = L.score.get(l.size);
  cl.stage(L.in.get(up), pl, crow_start);
  to_device(d, L.in.p, up, st);
  const Crops c = cl.at(d);
  label_chain(st, c, crows, cl.crop, words, cl.end - cl.nc, t);
  i32 *cnt = ptr<i32>(d, o_cnt), *rc = ptr<i32>(d, o_rc), *rh = ptr<i32>(d, o_rh), *kc = ptr<i32>(d, o_kc),
      *kh = ptr<i32>(d, o_kh);
  t.stage("topology", st, [&] {  // areas, hole counts and ranks in raster order
    zero(c.aux, sizeof(i32) * size_t(cl.crop), st);
    for_each(st, crows, Count{c, 1});
    for_each(st, crows, Count{c, 0});
    zero(cnt, sizeof(i32) * size_t(cl.crop), st);
    for_each(st, crows, Holes{c, cnt});
    zero(rc + crows, sizeof(i32), st), zero(rh + crows, sizeof(i32), st);  // each scan's trailing zero
    for_each(st, crows, RootCount{c, rc, rh});
    scan(st, rc, kc, crows + 1, L.scan);
    scan(st, rh, kh, crows + 1, L.scan);
    for_each(st, crows, CompEmit{c, kc, cnt, ptr<i32>(d, o_ca), ptr<i32>(d, o_ch)});
    for_each(st, crows, HoleEmit{c, kh, cnt, ptr<i32>(d, o_ha), ptr<i32>(d, o_hc)});
  });
  // The ranks of each plane's first crop row, then the per-component and per-hole numbers (and the crops).
  const size_t kb = sizeof(i32) * size_t(crows + 1), kb_pad = (kb + 255) / 256 * 256;
  u8* hk = L.out.get(kb_pad + kb);
  to_host(hk, kc, kb, st);
  to_host(hk + kb_pad, kh, kb, st);
  wait(st);
  const std::vector<i32> rk(ptr<i32>(hk, 0), ptr<i32>(hk, 0) + crows + 1),
      rhk(ptr<i32>(hk + kb_pad, 0), ptr<i32>(hk + kb_pad, 0) + crows + 1);
  const i64 C = rk.back(), Hn = rhk.back();
  Layout o;
  const size_t g_ca = o.put<i32>(C), g_ch = o.put<i32>(C), g_ha = o.put<i32>(Hn), g_hc = o.put<i32>(Hn),
               g_img = o.put<u8>(want_masks ? cl.crop : 0);
  u8* ho = L.out.get(o.size);
  to_host(ho + g_ca, d + o_ca, sizeof(i32) * size_t(C), st);
  to_host(ho + g_ch, d + o_ch, sizeof(i32) * size_t(C), st);
  to_host(ho + g_ha, d + o_ha, sizeof(i32) * size_t(Hn), st);
  to_host(ho + g_hc, d + o_hc, sizeof(i32) * size_t(Hn), st);
  if (want_masks) to_host(ho + g_img, c.img, size_t(cl.crop), st);
  wait(st);
  const i32 *ca = ptr<i32>(ho, g_ca), *ch = ptr<i32>(ho, g_ch), *ha = ptr<i32>(ho, g_ha), *hc = ptr<i32>(ho, g_hc);
  out.comp_area.assign(ca, ca + C);
  out.comp_holes.assign(ch, ch + C);
  for (i64 k = 0; k < C; ++k) out.hole_start.push_back(out.hole_start.back() + ch[k]);
  out.hole_area.resize(size_t(Hn));
  std::vector<i64> at(out.hole_start.begin(), out.hole_start.end() - 1);
  for (i64 j = 0; j < Hn; ++j) out.hole_area[size_t(at[size_t(hc[j])]++)] = ha[j];  // per component, raster order
  for (i64 k = 0; k < n; ++k) {
    const i64 a = rk[size_t(crow_start[size_t(k)])], b = rk[size_t(crow_start[size_t(k + 1)])];
    out.comp_start.push_back(b);
    std::vector<i32> s(ch + a, ch + b);
    std::sort(s.begin(), s.end(), [](i32 x, i32 y) { return x > y; });
    out.sig.insert(out.sig.end(), s.begin(), s.end());
    out.sig_start.push_back(i64(out.sig.size()));
    const Plane& p = pl[size_t(k)];
    if (!want_masks || !p.area) continue;
    const u8* img = ptr<u8>(ho, g_img) + p.crop;
    for (i32 y = 0; y < p.bh; ++y)
      std::copy_n(img + i64(y + 1) * (p.bw + 2) + 1, p.bw,
                  out.masks.data() + in.start[k] + i64(p.by + y) * p.W + p.bx);
  }
  return out;
}

}  // namespace gpucheck
