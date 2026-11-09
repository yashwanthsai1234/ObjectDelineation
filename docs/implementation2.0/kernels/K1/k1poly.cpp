// k1poly.cpp -- K1 polygonizer, single source in the style of REPO/cuda/gpucheck.cpp: every step is a __host__
// __device__ functor run by for_each over independent items (words, rows, runs, rings, seeds). for_each is a plain loop
// when the C++ compiler builds this file and a kernel launch when nvcc does (nvcc -x cu). Scans are a loop on the CPU
// and cub::DeviceScan::ExclusiveSum on the GPU. All arithmetic is integer, atomics are integer add / min / or only, so
// the result does not depend on thread order: the CPU build and the GPU build are bit-identical.
//
// Conventions (as gpucheck.cpp masks section and polygonize.py): vertex (x, y), 0 <= x <= W, 0 <= y <= H; edge = (vertex,
// direction d) with d 0 east, 1 south, 2 west, 3 north; foreground lies right of every edge; the bit-quad at vertex
// (x, y) is TL = pixel (x-1, y-1) = 1, TR = (x, y-1) = 2, BL = (x-1, y) = 4, BR = (x, y) = 8; at a saddle (codes 6, 9)
// the walk turns left: 8-connected foreground, 4-connected background.
//
// Steps of one chunk of masks (DESIGN.md has the proofs):
//   Pack      bytes -> bit words (pixel x = bit x % 32 of word x / 32, LSB first), one item per word.
//   Quads     bit-quad counts n1, n2, nD, n3, n4 by bit slicing on words, one item per vertex row.
//   RunCount  foreground / background runs per row from the words (popcount of run starts); scan -> run offsets.
//   [host sync 1: the run totals size every later array]
//   RunEmit   runs (x0, x1) per row, in raster order; background run 0 of every mask is the OUTSIDE sentinel.
//   Merge     union-find on runs: foreground rows y-1, y joined if 8-adjacent, background if 4-adjacent; background
//             runs on the image border joined to the sentinel. Roots are the minimum index = first run in raster order.
//   Flatten   every run points at its root.
//   HoleKey   per bounded background component B (= one hole ring) the key of its ring's first east edge in raster
//             order: min over fg pixels p of the surrounding component S(B) with the pixel above p in B.
//   RingMark  order bitmap (one bit per ring at its order pixel) and start bitmap (one bit at its start pixel).
//   scan      popcounts of both bitmaps -> ring index (= polygonize's order) and start rank, no sort.
//   RingInit  ring g's start edge in closed form: outer = (x0, y, east) of the component's first run, hole = (x0, y,
//             south) of the hole's first run.
//   Seeds     ring starts, plus (seed_rows > 0) the two vertical sides of every foreground run on rows y % K == 0.
//   Trace 1   each seed walks boundary edges (2x2 successor rule) to the next seed: corners, area, next seed.
//   RingWalk  per ring, the seed list in ring order: offset of each segment's first corner.
//   scan      ring corner counts -> output offsets.
//   Trace 2   each seed walks its segment again and writes its corners in place.
//   [host sync 2: ring and corner totals; sync 3: the outputs]
#include "k1poly.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <string>

#ifdef __CUDACC__
#include <cuda_runtime.h>

#include <cub/device/device_scan.cuh>
#define HD __host__ __device__
#else
#define HD
#endif

namespace k1poly {
namespace {

// ===================================================================== platform
#ifdef __CUDACC__
void cuda_check(cudaError_t e) {
  if (e != cudaSuccess) throw std::runtime_error(std::string("CUDA: ") + cudaGetErrorString(e));
}
using Stream = cudaStream_t;
Stream new_stream() {
  cudaStream_t s = nullptr;
  cuda_check(cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking));
  return s;
}
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
__global__ void kernel(long long n, F f) {
  const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) f(i64(i));
}
template <class F>
void for_each(Stream st, i64 n, const F& f) {
  constexpr int kThreads = 128;
  if (n <= 0) return;
  kernel<<<unsigned((n + kThreads - 1) / kThreads), kThreads, 0, st>>>((long long)n, f);
  cuda_check(cudaGetLastError());
}
#else
struct Stream {};
Stream new_stream() { return {}; }
void* alloc(bool, size_t bytes) { return ::operator new(bytes); }
void release(bool, void* p) { ::operator delete(p); }
void to_device(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void to_host(void* dst, const void* src, size_t bytes, Stream) { if (bytes) std::memcpy(dst, src, bytes); }
void zero(void* p, size_t bytes, Stream) { if (bytes) std::memset(p, 0, bytes); }
void wait(Stream) {}
// K1_ORDER (tests only) runs the items in another order to stand in for the GPU's arbitrary schedule: 0 forward,
// 1 reverse, 2 the permutation k -> k p mod n (p = 2^31 - 1, prime, so coprime with any n < p). Every step must give
// the same result in every order.
#ifndef K1_ORDER
#define K1_ORDER 0
#endif
template <class F>
void for_each(Stream, i64 n, const F& f) {
  for (i64 k = 0; k < n; ++k) {
    const i64 i = K1_ORDER == 0 ? k : K1_ORDER == 1 ? n - 1 - k : (k * 2147483647LL) % n;
    f(i);
  }
}
#endif

template <bool kDevice>
struct Block {  // grow-only buffer (next power of two), never freed per call
  unsigned char* p = nullptr;
  size_t cap = 0;
  unsigned char* get(size_t bytes) {
    if (bytes <= cap) return p;
    size_t c = 4096;
    while (c < bytes) c <<= 1;
    if (p) release(kDevice, p);
    p = nullptr, cap = 0;
    p = static_cast<unsigned char*>(alloc(kDevice, c));
    cap = c;
    return p;
  }
};

// Exclusive prefix sum of n items (out[k] = in[0] + ... + in[k-1]); callers pass one trailing zero to get the total.
void scan(Stream st, const i32* in, i32* out, i64 n, Block<true>& tmp) {
  if (n <= 0) return;
#ifdef __CUDACC__
  size_t bytes = 0;
  cuda_check(cub::DeviceScan::ExclusiveSum(nullptr, bytes, in, out, int(n), st));
  void* t = tmp.get(bytes > 0 ? bytes : 1);
  cuda_check(cub::DeviceScan::ExclusiveSum(t, bytes, in, out, int(n), st));
#else
  (void)st, (void)tmp;
  i32 s = 0;
  for (i64 k = 0; k < n; ++k) {
    const i32 v = in[k];
    out[k] = s;
    s += v;
  }
#endif
}

struct Engine {
  std::mutex mu;
  bool ready = false;
  Stream stream{};
  Block<true> devA, devB, tmp;
  Block<false> in, out;
  void start() {
    if (!ready) stream = new_stream(), ready = true;
  }
};
Engine& engine() {
  static Engine* e = new Engine;  // never destroyed (the CUDA runtime may be torn down first)
  return *e;
}

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
T* ptr(unsigned char* base, size_t off) { return reinterpret_cast<T*>(base + off); }

// Stage clock: CUDA events on the stream (read after the last sync), wall clock on the CPU build.
struct Clock {
#ifdef __CUDACC__
  cudaEvent_t ev[8];
  Clock() { for (auto& e : ev) cuda_check(cudaEventCreate(&e)); }
  ~Clock() { for (auto& e : ev) cudaEventDestroy(e); }
  void mark(int k, Stream st) { cuda_check(cudaEventRecord(ev[k], st)); }
  double ms(int a, int b) { float t = 0; cuda_check(cudaEventElapsedTime(&t, ev[a], ev[b])); return t; }
#else
  std::chrono::steady_clock::time_point t[8];
  void mark(int k, Stream) { t[k] = std::chrono::steady_clock::now(); }
  double ms(int a, int b) { return std::chrono::duration<double, std::milli>(t[b] - t[a]).count(); }
#endif
};

// ===================================================================== atomics and bit helpers
HD inline void atomic_add(i64* p, i64 v) {
#ifdef __CUDA_ARCH__
  atomicAdd(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(v));
#else
  *p += v;
#endif
}
HD inline void atomic_min(i32* p, i32 v) {
#ifdef __CUDA_ARCH__
  atomicMin(p, v);
#else
  if (v < *p) *p = v;
#endif
}
HD inline i32 atomic_min_old(i32* p, i32 v) {
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
HD inline void set_flag(i32* p) {  // error flag: any thread may set it
#ifdef __CUDA_ARCH__
  atomicExch(p, 1);
#else
  *p = 1;
#endif
}
HD inline i32 load(const i32* p) {  // a fresh read of union-find parents other threads write
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
  return __ffs(int(x)) - 1;
#else
  return __builtin_ctz(x);
#endif
}

// Segment of item i in a prefix array start[0..n]: the largest s < n with start[s] <= i (empty segments skipped).
HD inline i32 segment(const i32* start, i32 n, i64 i) {
  i32 lo = 0, hi = n;
  while (hi - lo > 1) {
    const i32 mid = (lo + hi) / 2;
    if (start[mid] <= i) lo = mid;
    else hi = mid;
  }
  return lo;
}

HD inline i32 step_x(int d) { return d == 0 ? 1 : (d == 2 ? -1 : 0); }
HD inline i32 step_y(int d) { return d == 1 ? 1 : (d == 3 ? -1 : 0); }
HD inline bool has_edge(int code, int d) {  // gpucheck.cpp has_edge
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

// ===================================================================== device views
struct Geo {  // one chunk of masks
  const i32 *W, *H, *wpr;  // per mask: width, height, words per row
  const i32 *ws, *rs, *vs, *ss;  // prefix (n + 1): words, rows, vertex rows (H + 1), seed-map entries
  const i32* rowmask;      // per row: its mask
  i32 n, K;                // masks, seed rows (0 = off)

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
  // straight, else right (gpucheck.cpp InitSlots, polygonize.link).
  HD int succ_dir(const u32* bits, i32 m, i32 x, i32 y, int d) const {
    const int c = quad(bits, m, x, y);
    if (has_edge(c, (d + 3) & 3)) return (d + 3) & 3;
    if (has_edge(c, d)) return d;
    if (has_edge(c, (d + 1) & 3)) return (d + 1) & 3;
    return -1;
  }
  // Direction of the edge arriving at vertex (x, y) whose successor leaves in direction t (the successor map is a
  // permutation of the edges, so exactly one exists).
  HD int pred_dir(const u32* bits, i32 m, i32 x, i32 y, int t) const {
    for (int d = 0; d < 4; ++d) {
      const i32 ux = x - step_x(d), uy = y - step_y(d);
      if (has_edge(quad(bits, m, ux, uy), d) && succ_dir(bits, m, x, y, d) == t) return d;
    }
    return -1;
  }
  HD u32 valid(i32 m, i32 j) const {  // bits of word j inside the row
    const i32 r = W[m] - 32 * j;
    return r >= 32 ? ~0u : (r <= 0 ? 0u : ((1u << r) - 1u));
  }
  HD i32 bitword(i32 m, i32 x, i32 y) const { return ws[m] + y * wpr[m] + (x >> 5); }
};

struct Runs {
  const i32 *rf, *rb;            // per row (rows + 1): first foreground / background run
  i32 *fx0, *fx1, *frow, *fpar;  // foreground runs
  i32 *bx0, *bx1, *brow, *bpar, *bleft, *hkey;  // background runs; bx0 = -1 marks a mask's OUTSIDE sentinel
};

HD inline i32 find(const i32* par, i32 x) {
  for (;;) {
    const i32 p = load(par + x);
    if (p == x) return x;
    x = p;
  }
}

// Lock-free union with the root = minimum index (Playne and Hawick 2018, Alg. "merge"; also the union of
// Allegretti et al. 2019 and gpucheck.cpp's scorer labeling).
HD inline void unite(i32* par, i32 a, i32 b) {
  for (;;) {
    a = find(par, a);
    b = find(par, b);
    if (a == b) return;
    if (a < b) {
      const i32 old = atomic_min_old(par + b, a);
      if (old == b) return;
      b = old;
    } else {
      const i32 old = atomic_min_old(par + a, b);
      if (old == a) return;
      a = old;
    }
  }
}

// ===================================================================== stage 1: words, quads, runs
struct Pack {
  Geo g;
  const u8* bytes;
  const i64* start;  // chunk-local byte offset of each mask
  u32* words;
  HD void operator()(i64 w) const {
    const i32 m = segment(g.ws, g.n, w);
    const i32 lw = i32(w - g.ws[m]), wpr = g.wpr[m], W = g.W[m];
    const i32 y = lw / wpr, j = lw - y * wpr, x0 = 32 * j;
    const u8* row = bytes + start[m] + i64(y) * W;
    u32 v = 0;
    const i32 nb = W - x0 < 32 ? W - x0 : 32;
    for (i32 b = 0; b < nb; ++b) v |= u32(row[x0 + b] != 0) << b;
    words[w] = v;
  }
};

#ifdef __CUDACC__
// GPU packing: one warp per word, lane l reads pixel 32 j + l (coalesced bytes) and __ballot_sync builds the word
// (as gpucheck.cpp's raster). Same words as Pack. The early exit is uniform per warp (blockDim is a multiple of 32
// and all 32 lanes of a warp share w), so every lane named in the ballot mask executes it.
__global__ void pack_kernel(Pack p, long long nw) {
  const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x, w = t >> 5;
  const int lane = int(threadIdx.x & 31u);
  if (w >= nw) return;
  const Geo& g = p.g;
  const i32 m = segment(g.ws, g.n, w);
  const i32 lw = i32(w - g.ws[m]), wpr = g.wpr[m], W = g.W[m];
  const i32 y = lw / wpr, x = 32 * (lw - y * wpr) + lane;
  const bool v = x < W && p.bytes[p.start[m] + i64(y) * W + x] != 0;
  const u32 word = __ballot_sync(0xffffffffu, v);
  if (lane == 0) p.words[w] = word;
}
void pack(Stream st, const Pack& p, i64 nw) {
  if (nw <= 0) return;
  constexpr int kThreads = 256;
  const long long threads = (long long)nw * 32;
  pack_kernel<<<unsigned((threads + kThreads - 1) / kThreads), kThreads, 0, st>>>(p, (long long)nw);
  cuda_check(cudaGetLastError());
}
#else
void pack(Stream st, const Pack& p, i64 nw) { for_each(st, nw, p); }
#endif

// Bit-quad counts by bit slicing over 32 vertex columns at once (checker.bitquad: n1, n2, nD, n3, n4).
struct Quads {
  Geo g;
  const u32* words;
  i64* quads;
  HD void operator()(i64 i) const {
    const i32 m = segment(g.vs, g.n, i);
    const i32 y = i32(i - g.vs[m]), W = g.W[m];
    const i32 nvw = (W + 1 + 31) / 32;
    i64 n1 = 0, n2 = 0, nd = 0, n3 = 0, n4 = 0;
    for (i32 j = 0; j < nvw; ++j) {
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
    i32 fi = R.rf[r], bi = R.rb[r];
    if (y == 0) {
      R.bx0[bi] = -1, R.bx1[bi] = -1, R.brow[bi] = i32(r), R.bpar[bi] = bi, R.bleft[bi] = -1, R.hkey[bi] = INT_MAX;
      ++bi;
    }
    const i32 f0 = fi, b0 = bi;
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
    {  // foreground, 8-adjacent: a.x0 <= b.x1 + 1 and b.x0 <= a.x1 + 1
      i32 a = R.rf[r], ae = R.rf[r + 1], b = R.rf[r - 1], bend = R.rf[r];
      while (a < ae && b < bend) {
        if (R.fx0[a] <= R.fx1[b] + 1 && R.fx0[b] <= R.fx1[a] + 1) unite(R.fpar, a, b);
        const i32 xa = R.fx1[a], xb = R.fx1[b];
        if (xa <= xb) ++a;
        if (xb <= xa) ++b;
      }
    }
    {  // background, 4-adjacent: a.x0 <= b.x1 and b.x0 <= a.x1
      i32 a = b_lo, ae = b_hi, b = R.rb[r - 1] + (y - 1 == 0), bend = R.rb[r];
      while (a < ae && b < bend) {
        if (R.bx0[a] <= R.bx1[b] && R.bx0[b] <= R.bx1[a]) unite(R.bpar, a, b);
        const i32 xa = R.bx1[a], xb = R.bx1[b];
        if (xa <= xb) ++a;
        if (xb <= xa) ++b;
      }
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

// Order key of each hole ring: its first east edge in raster order. An east edge of the hole ring of B is the top
// side of a foreground pixel p of S(B) (the component surrounding B) whose upper neighbour is in B. S(B) is the
// component of the foreground run just left of B's first run (DESIGN.md, lemma 3).
struct HoleKey {
  Geo g;
  Runs R;
  i32* err;
  HD void operator()(i64 r) const {
    const i32 m = g.rowmask[r], y = i32(r - g.rs[m]), W = g.W[m];
    if (y == 0) return;
    i32 a = R.rf[r], ae = R.rf[r + 1], b = R.rb[r - 1] + (y - 1 == 0), bend = R.rb[r];
    while (a < ae && b < bend) {
      if (R.fx0[a] <= R.bx1[b] && R.bx0[b] <= R.fx1[a]) {
        const i32 B = R.bpar[b];
        if (R.bx0[B] >= 0) {  // a hole (not the outside)
          const i32 left = R.bleft[B];
          if (left < 0) set_flag(err);
          else if (R.fpar[a] == R.fpar[left]) {
            const i32 x = R.fx0[a] > R.bx0[b] ? R.fx0[a] : R.bx0[b];
            atomic_min(R.hkey + B, y * W + x);
          }
        }
      }
      const i32 xa = R.fx1[a], xb = R.bx1[b];
      if (xa <= xb) ++a;
      if (xb <= xa) ++b;
    }
  }
};

// One bit per ring in the order bitmap (at its order pixel) and in the start bitmap (at its start pixel).
struct RingMark {
  Geo g;
  Runs R;
  i32 nf;
  u32 *obits, *sbits;
  i32* err;
  HD void operator()(i64 i) const {
    if (i < nf) {
      if (R.fpar[i] != i) return;
      const i32 r = R.frow[i], m = g.rowmask[r], y = r - g.rs[m], x = R.fx0[i];
      atomic_or(obits + g.bitword(m, x, y), 1u << (x & 31));
      atomic_or(sbits + g.bitword(m, x, y), 1u << (x & 31));
    } else {
      const i32 b = i32(i - nf);
      if (R.bpar[b] != b || R.bx0[b] < 0) return;
      const i32 r = R.brow[b], m = g.rowmask[r], y = r - g.rs[m], x = R.bx0[b], W = g.W[m];
      atomic_or(sbits + g.bitword(m, x, y), 1u << (x & 31));
      const i32 k = R.hkey[b];
      if (k == INT_MAX) { set_flag(err); return; }
      const i32 ky = k / W, kx = k - ky * W;
      atomic_or(obits + g.bitword(m, kx, ky), 1u << (kx & 31));
    }
  }
};

struct Popc {
  const u32* bits;
  i32* cnt;
  HD void operator()(i64 w) const { cnt[w] = popc(bits[w]); }
};

HD inline i32 rank_of(const u32* bits, const i32* scan, i32 w, i32 x) {
  return scan[w] + popc(bits[w] & ((1u << (x & 31)) - 1u));
}

struct Seeds {
  i32 *x, *y, *d, *m;            // start edge of each seed; seeds [0, R) are the ring starts, in ring order
  i32 *cnt, *next, *off, *ring, *dead;
  i64* area;
};

struct RingInit {
  Geo g;
  Runs R;
  i32 nf;
  const u32 *obits, *sbits;
  const i32 *oscan, *sscan;
  i32* ring_of_start;
  Seeds S;
  HD void operator()(i64 i) const {
    i32 r, x, ox, oy, d;
    if (i < nf) {
      if (R.fpar[i] != i) return;
      r = R.frow[i], x = R.fx0[i], d = 0;
    } else {
      const i32 b = i32(i - nf);
      if (R.bpar[b] != b || R.bx0[b] < 0) return;
      r = R.brow[b], x = R.bx0[b], d = 1;
    }
    const i32 m = g.rowmask[r], y = r - g.rs[m];
    if (d == 0) ox = x, oy = y;
    else oy = R.hkey[i - nf] / g.W[m], ox = R.hkey[i - nf] - oy * g.W[m];
    const i32 ring = rank_of(obits, oscan, g.bitword(m, ox, oy), ox);
    ring_of_start[rank_of(sbits, sscan, g.bitword(m, x, y), x)] = ring;
    S.x[ring] = x, S.y[ring] = y, S.d[ring] = d, S.m[ring] = m;
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
    S.x[id] = x0, S.y[id] = y + 1, S.d[id] = 3, S.m[id] = m;          // left side of pixel (x0, y), going north
    S.x[id + 1] = x1 + 1, S.y[id + 1] = y, S.d[id + 1] = 1, S.m[id + 1] = m;  // right side of pixel (x1, y), south
  }
};

struct SetCounts {  // counts[0] = rings R, [1] = seeds; Header sets [2] = corners; [4] = error flag
  const i32 *oscan, *soff;
  i32 nw, nf, K;
  i32* counts;
  HD void operator()(i64) const {
    counts[0] = oscan[nw];
    counts[1] = counts[0] + (K > 0 ? 2 * soff[nf] : 0);
  }
};

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
    const i32 s = i32(i), R = counts[0];
    if (s >= counts[1]) return;
    const i32 m = S.m[s];
    i32 x = S.x[s], y = S.y[s];
    int d = S.d[s];
    if (pass == 1) {
      if (s >= R && seed_of(m, x, y, d) != s) {  // a row seed on a ring start edge: the ring start covers it
        S.dead[s] = 1;
        return;
      }
    } else if (S.dead[s]) {
      return;
    }
    int prev = g.pred_dir(words, m, x, y, d);
    i64 k = pass == 2 ? i64(rbase[S.ring[s]]) + S.off[s] : 0;
    i32 corners = 0, nxt = -1;
    i64 area = 0;
    const i64 cap = 4 * (i64(g.W[m]) + 1) * (i64(g.H[m]) + 1);
    for (i64 step = 0;; ++step) {
      if (step > cap || d < 0) { set_flag(err); return; }
      if (d != prev) {
        if (pass == 2) pts[2 * k] = x, pts[2 * k + 1] = y, ++k;
        ++corners;
      }
      area += area_term(x, y, d);
      x += step_x(d), y += step_y(d);
      prev = d;
      d = g.succ_dir(words, m, x, y, d);
      if (d < 0) { set_flag(err); return; }
      nxt = seed_of(m, x, y, d);
      if (nxt >= 0) break;
    }
    if (pass == 1) S.cnt[s] = corners, S.area[s] = area, S.next[s] = nxt;
  }
};

struct RingWalk {
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
      if (++n > cap || s < 0) { set_flag(err); return; }
    } while (s != g0);
    rcount[g0] = off, rarea[g0] = a;
  }
};

struct Header {  // counts[2] = corners total; ring_start of each mask from the order scan
  Geo g;
  const i32 *oscan, *rbase;
  i32 rmax;
  i32 *counts, *ring_start;
  HD void operator()(i64 m) const {
    if (m == 0) counts[2] = rbase[rmax];
    ring_start[m] = oscan[g.ws[m]];
  }
};

// ===================================================================== one chunk
// Returns false (nothing appended) when the chunk's run-sized scratch would exceed opt.max_scratch_bytes and the chunk
// has more than one mask: the caller then splits it in halves (the steps before sync 1 are cheap to redo).
bool chunk(const Masks& in, i64 m0, i64 m1, const Options& opt, Result& out, Engine& E, Timing* tm) {
  const i32 n = i32(m1 - m0), K = opt.seed_rows;
  std::vector<i32> W(n), H(n), wpr(n), ws(n + 1, 0), rs(n + 1, 0), vs(n + 1, 0), ss(n + 1, 0);
  std::vector<i64> start(n + 1, 0);
  for (i32 k = 0; k < n; ++k) {
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

  // ---- layout A: inputs, words, bitmaps, row arrays
  Layout la;
  const size_t o_bytes = la.put<u8>(npx), o_start = la.put<i64>(n + 1), o_W = la.put<i32>(n), o_H = la.put<i32>(n),
               o_wpr = la.put<i32>(n), o_ws = la.put<i32>(n + 1), o_rs = la.put<i32>(n + 1), o_vs = la.put<i32>(n + 1),
               o_ss = la.put<i32>(n + 1), o_rowmask = la.put<i32>(NR);
  const size_t in_bytes = la.size;
  const size_t o_words = la.put<u32>(NW), o_cf = la.put<i32>(NR + 1), o_cb = la.put<i32>(NR + 1),
               o_rf = la.put<i32>(NR + 1), o_rb = la.put<i32>(NR + 1), o_ocnt = la.put<i32>(NW + 1),
               o_scnt = la.put<i32>(NW + 1), o_oscan = la.put<i32>(NW + 1), o_sscan = la.put<i32>(NW + 1),
               o_seedmap = la.put<i32>(ss[n]);
  const size_t z0 = la.size;  // [z0, end) zeroed
  const size_t o_obits = la.put<u32>(NW), o_sbits = la.put<u32>(NW), o_quads = la.put<i64>(5 * i64(n)),
               o_counts = la.put<i32>(8);
  unsigned char* A = E.devA.get(la.size);
  unsigned char* h = E.in.get(std::max(in_bytes, size_t(64)));
  if (npx) std::memcpy(ptr<u8>(h, o_bytes), in.bits + in.start[m0], size_t(npx));
  std::memcpy(ptr<i64>(h, o_start), start.data(), sizeof(i64) * (n + 1));
  std::memcpy(ptr<i32>(h, o_W), W.data(), sizeof(i32) * n);
  std::memcpy(ptr<i32>(h, o_H), H.data(), sizeof(i32) * n);
  std::memcpy(ptr<i32>(h, o_wpr), wpr.data(), sizeof(i32) * n);
  std::memcpy(ptr<i32>(h, o_ws), ws.data(), sizeof(i32) * (n + 1));
  std::memcpy(ptr<i32>(h, o_rs), rs.data(), sizeof(i32) * (n + 1));
  std::memcpy(ptr<i32>(h, o_vs), vs.data(), sizeof(i32) * (n + 1));
  std::memcpy(ptr<i32>(h, o_ss), ss.data(), sizeof(i32) * (n + 1));
  if (NR) std::memcpy(ptr<i32>(h, o_rowmask), rowmask.data(), sizeof(i32) * NR);

  const Stream st = E.stream;
  Clock clk;
  clk.mark(0, st);
  to_device(A, h, in_bytes, st);
  zero(A + z0, la.size - z0, st);
  zero(A + o_cf, sizeof(i32) * (NR + 1), st);
  zero(A + o_cb, sizeof(i32) * (NR + 1), st);
  zero(A + o_ocnt, sizeof(i32) * (NW + 1), st);
  zero(A + o_scnt, sizeof(i32) * (NW + 1), st);
  clk.mark(1, st);
  const Geo g{ptr<i32>(A, o_W), ptr<i32>(A, o_H), ptr<i32>(A, o_wpr), ptr<i32>(A, o_ws), ptr<i32>(A, o_rs),
              ptr<i32>(A, o_vs), ptr<i32>(A, o_ss), ptr<i32>(A, o_rowmask), n, K};
  u32* words = ptr<u32>(A, o_words);
  i32* counts = ptr<i32>(A, o_counts);
  i32* err = counts + 4;
  pack(st, Pack{g, ptr<u8>(A, o_bytes), ptr<i64>(A, o_start), words}, NW);
  for_each(st, NV, Quads{g, words, ptr<i64>(A, o_quads)});
  for_each(st, NR, RunCount{g, words, ptr<i32>(A, o_cf), ptr<i32>(A, o_cb)});
  scan(st, ptr<i32>(A, o_cf), ptr<i32>(A, o_rf), NR + 1, E.tmp);
  scan(st, ptr<i32>(A, o_cb), ptr<i32>(A, o_rb), NR + 1, E.tmp);
  unsigned char* o = E.out.get(64);
  to_host(o, A + o_rf + sizeof(i32) * NR, sizeof(i32), st);
  to_host(o + 4, A + o_rb + sizeof(i32) * NR, sizeof(i32), st);
  wait(st);  // sync 1
  const i32 NF = ptr<i32>(o, 0)[0], NB = ptr<i32>(o, 4)[0];
  const i64 rmax64 = i64(NF) + NB, smax64 = rmax64 + 2 * i64(NF);
  if (NF < 0 || NB < 0 || smax64 > INT_MAX / 2 || 8 * i64(NF) > INT_MAX)
    throw std::length_error("k1poly: chunk too large (lower Options::chunk_pixels)");
  const i32 RMAX = i32(rmax64), SMAX = i32(smax64);

  // ---- layout B: runs, rings, seeds, corners
  Layout lb;
  size_t o_f[4], o_b[6];
  for (size_t& q : o_f) q = lb.put<i32>(NF);
  for (size_t& q : o_b) q = lb.put<i32>(NB);
  const size_t o_flag = lb.put<i32>(NF + 1), o_soff = lb.put<i32>(NF + 1), o_ros = lb.put<i32>(RMAX),
               o_rbase = lb.put<i32>(RMAX + 1), o_rarea = lb.put<i64>(RMAX), o_pts = lb.put<i32>(8 * i64(NF)),
               o_rstart = lb.put<i32>(n + 1);
  size_t o_s[8];
  for (size_t& q : o_s) q = lb.put<i32>(SMAX);
  const size_t o_sarea = lb.put<i64>(SMAX);
  const size_t zb = lb.size;  // [zb, end) zeroed
  const size_t o_rcount = lb.put<i32>(RMAX + 1), o_dead = lb.put<i32>(SMAX);
  if (n > 1 && i64(lb.size) > opt.max_scratch_bytes) return false;
  unsigned char* B = E.devB.get(lb.size);
  zero(B + zb, lb.size - zb, st);
  zero(B + o_flag, sizeof(i32) * (NF + 1), st);
  const Runs R{ptr<i32>(A, o_rf), ptr<i32>(A, o_rb), ptr<i32>(B, o_f[0]), ptr<i32>(B, o_f[1]), ptr<i32>(B, o_f[2]),
               ptr<i32>(B, o_f[3]), ptr<i32>(B, o_b[0]), ptr<i32>(B, o_b[1]), ptr<i32>(B, o_b[2]), ptr<i32>(B, o_b[3]),
               ptr<i32>(B, o_b[4]), ptr<i32>(B, o_b[5])};
  const Seeds S{ptr<i32>(B, o_s[0]), ptr<i32>(B, o_s[1]), ptr<i32>(B, o_s[2]), ptr<i32>(B, o_s[3]),
                ptr<i32>(B, o_s[4]), ptr<i32>(B, o_s[5]), ptr<i32>(B, o_s[6]), ptr<i32>(B, o_s[7]),
                ptr<i32>(B, o_dead), ptr<i64>(B, o_sarea)};
  u32 *obits = ptr<u32>(A, o_obits), *sbits = ptr<u32>(A, o_sbits);
  i32 *oscan = ptr<i32>(A, o_oscan), *sscan = ptr<i32>(A, o_sscan), *ros = ptr<i32>(B, o_ros);
  i32 *rcount = ptr<i32>(B, o_rcount), *rbase = ptr<i32>(B, o_rbase);

  for_each(st, NR, RunEmit{g, words, R});
  for_each(st, NR, Merge{g, R});
  for_each(st, rmax64, Flatten{R, NF});
  clk.mark(2, st);
  for_each(st, NR, HoleKey{g, R, err});
  for_each(st, rmax64, RingMark{g, R, NF, obits, sbits, err});
  for_each(st, NW, Popc{obits, ptr<i32>(A, o_ocnt)});
  for_each(st, NW, Popc{sbits, ptr<i32>(A, o_scnt)});
  scan(st, ptr<i32>(A, o_ocnt), oscan, NW + 1, E.tmp);
  scan(st, ptr<i32>(A, o_scnt), sscan, NW + 1, E.tmp);
  for_each(st, rmax64, RingInit{g, R, NF, obits, sbits, oscan, sscan, ros, S});
  if (K > 0) {
    for_each(st, NF, SeedFlag{g, R, ptr<i32>(B, o_flag)});
    scan(st, ptr<i32>(B, o_flag), ptr<i32>(B, o_soff), NF + 1, E.tmp);
  }
  for_each(st, 1, SetCounts{oscan, ptr<i32>(B, o_soff), NW, NF, K, counts});
  if (K > 0) for_each(st, NF, SeedEmit{g, R, ptr<i32>(B, o_flag), ptr<i32>(B, o_soff), counts, ptr<i32>(A, o_seedmap), S});
  clk.mark(3, st);
  Trace t1{g, words, sbits, sscan, ros, ptr<i32>(A, o_seedmap), counts, rbase, S, ptr<i32>(B, o_pts), err, 1};
  for_each(st, SMAX, t1);
  for_each(st, RMAX, RingWalk{counts, S, rcount, ptr<i64>(B, o_rarea), err, SMAX});
  scan(st, rcount, rbase, RMAX + 1, E.tmp);
  Trace t2 = t1;
  t2.pass = 2;
  for_each(st, SMAX, t2);
  for_each(st, n + 1, Header{g, oscan, rbase, RMAX, counts, ptr<i32>(B, o_rstart)});
  clk.mark(4, st);
  to_host(o, counts, sizeof(i32) * 8, st);
  wait(st);  // sync 2
  i32 hdr[8];
  std::memcpy(hdr, o, sizeof(hdr));
  if (hdr[4]) throw std::logic_error("k1poly: invariant violated on the device (trace, hole key or ring walk)");
  const i32 nr = hdr[0], ncorner = hdr[2];
  if (nr < 0 || nr > RMAX || ncorner < 0 || ncorner > 4 * NF) throw std::logic_error("k1poly: counts out of bound");
  Layout lo;
  const size_t q_rbase = lo.put<i32>(nr + 1), q_area = lo.put<i64>(nr), q_pts = lo.put<i32>(2 * i64(ncorner)),
               q_rstart = lo.put<i32>(n + 1), q_quads = lo.put<i64>(5 * i64(n));
  o = E.out.get(lo.size);
  to_host(o + q_rbase, rbase, sizeof(i32) * (nr + 1), st);
  to_host(o + q_area, ptr<i64>(B, o_rarea), sizeof(i64) * nr, st);
  to_host(o + q_pts, ptr<i32>(B, o_pts), sizeof(i32) * 2 * size_t(ncorner), st);
  to_host(o + q_rstart, ptr<i32>(B, o_rstart), sizeof(i32) * (n + 1), st);
  to_host(o + q_quads, ptr<i64>(A, o_quads), sizeof(i64) * 5 * n, st);
  clk.mark(5, st);
  wait(st);  // sync 3
  if (tm) {
    tm->h2d += clk.ms(0, 1), tm->ccl += clk.ms(1, 2), tm->rings += clk.ms(2, 3), tm->trace += clk.ms(3, 4);
    tm->d2h += clk.ms(4, 5), tm->total += clk.ms(0, 5);
  }
  const i64 r0 = i64(out.ring_area.size()), p0 = out.pts_start.back();
  const i32* rst = ptr<i32>(o, q_rstart);
  for (i32 k = 1; k <= n; ++k) out.ring_start.push_back(r0 + rst[k]);
  const i32* rb = ptr<i32>(o, q_rbase);
  const i64* ar = ptr<i64>(o, q_area);
  for (i32 k = 0; k < nr; ++k) out.pts_start.push_back(p0 + rb[k + 1]), out.ring_area.push_back(ar[k] / 2);
  const i32* pp = ptr<i32>(o, q_pts);
  out.pts.insert(out.pts.end(), pp, pp + 2 * size_t(ncorner));
  const i64* qq = ptr<i64>(o, q_quads);
  out.quads.insert(out.quads.end(), qq, qq + 5 * size_t(n));
  return true;
}

void run_range(const Masks& in, i64 m0, i64 m1, const Options& opt, Result& out, Engine& E, Timing* tm) {
  if (chunk(in, m0, m1, opt, out, E, tm)) return;
  const i64 mid = m0 + (m1 - m0) / 2;
  run_range(in, m0, mid, opt, out, E, tm);
  run_range(in, mid, m1, opt, out, E, tm);
}

}  // namespace

bool cuda_build() {
#ifdef __CUDACC__
  return true;
#else
  return false;
#endif
}

std::string device_info() {
#ifdef __CUDACC__
  int dev = 0;
  cudaDeviceProp p{};
  cuda_check(cudaGetDevice(&dev));
  cuda_check(cudaGetDeviceProperties(&p, dev));
  return std::string(p.name) + ", " + std::to_string(p.multiProcessorCount) + " SMs, L2 " +
         std::to_string(p.l2CacheSize >> 20) + " MiB, CC " + std::to_string(p.major) + "." + std::to_string(p.minor) +
         ", " + std::to_string(p.totalGlobalMem >> 20) + " MiB";
#else
  return "cpu";
#endif
}

Result polygonize_masks(const Masks& in, const Options& opt, Timing* timing) {
  Result out;
  out.ring_start = {0};
  out.pts_start = {0};
  if (opt.seed_rows < 0) throw std::invalid_argument("seed_rows must be >= 0");
  for (i64 m = 0; m < in.n; ++m)
    if (in.W[m] < 0 || in.H[m] < 0 || in.W[m] > kMaxDim || in.H[m] > kMaxDim ||
        in.start[m + 1] - in.start[m] != i64(in.W[m]) * in.H[m])
      throw std::invalid_argument("every mask must be (H, W) with 0 <= H, W <= 4096 and H * W pixels");
  if (in.n == 0) return out;
  Engine& E = engine();
  std::lock_guard<std::mutex> lock(E.mu);
  E.start();
  for (i64 m0 = 0; m0 < in.n;) {
    i64 m1 = m0, px = 0;
    while (m1 < in.n) {
      const i64 p = i64(in.W[m1]) * in.H[m1];
      if (m1 > m0 && px + p > opt.chunk_pixels) break;
      px += p, ++m1;
    }
    run_range(in, m0, m1, opt, out, E, timing);
    m0 = m1;
  }
  return out;
}

}  // namespace k1poly
