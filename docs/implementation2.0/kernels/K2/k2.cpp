// k2.cpp -- host driver of the K2 kernels. It includes gpucheck.cpp (the snapshot in vendor/, unmodified) and then
// k2score.inc inside gpucheck's unnamed namespace, so the K2 functors run beside the existing ones on the same
// platform layer. Build: the C++ compiler (CPU build, -ffp-contract=off) or nvcc -x cu (GPU build, --fmad=false).
// The drivers are simple on purpose (one stream, grow-only device arenas, pageable copies): the production path is to
// paste k2score.inc into gpucheck.cpp and call the functors from its own Lane / Layout machinery (DESIGN.md).
#include "vendor/gpucheck.cpp"
#include "k2.h"

#include <chrono>
#include <cmath>
#include <limits>

namespace gpucheck {
namespace {
#include "k2score.inc"

#ifdef __CUDACC__
struct Timer {
  Stream st;
  cudaEvent_t a{}, b{};
  explicit Timer(Stream s) : st(s) {
    cuda_check(cudaEventCreate(&a));
    cuda_check(cudaEventCreate(&b));
  }
  ~Timer() {
    cudaEventDestroy(a);
    cudaEventDestroy(b);
  }
  void start() { cuda_check(cudaEventRecord(a, st)); }
  double stop() {
    cuda_check(cudaEventRecord(b, st));
    cuda_check(cudaEventSynchronize(b));
    float ms = 0;
    cuda_check(cudaEventElapsedTime(&ms, a, b));
    return ms;
  }
};
#else
struct Timer {
  std::chrono::steady_clock::time_point t0;
  explicit Timer(Stream) {}
  void start() { t0 = std::chrono::steady_clock::now(); }
  double stop() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(); }
};
#endif

Stream k2_stream() {
  static Stream s = new_stream();
  return s;
}

template <class T>
void put(u8* dev, size_t off, const std::vector<T>& v, Stream st) {
  to_device(dev + off, v.data(), sizeof(T) * v.size(), st);
}

template <class T>
std::vector<T> get(const u8* dev, size_t off, i64 n, Stream st) {
  std::vector<T> v(size_t(n > 0 ? n : 0));
  to_host(v.data(), dev + off, sizeof(T) * v.size(), st);
  wait(st);
  return v;
}

// Masks of a Bits batch plus their derived starts.
struct HostFrames {
  i64 n = 0, words = 0, pix = 0, rows = 0, cols = 0;
  std::vector<i64> word_start, pix_start, row_start, col_start;
  std::vector<i32> H, W;
  explicit HostFrames(const k2::Bits& b) : n(b.n), word_start(b.word_start, b.word_start + b.n + 1),
                                          pix_start(size_t(b.n + 1), 0), row_start(size_t(b.n + 1), 0),
                                          col_start(size_t(b.n + 1), 0), H(b.H, b.H + b.n), W(b.W, b.W + b.n) {
    for (i64 m = 0; m < n; ++m) {
      if (word_start[size_t(m + 1)] - word_start[size_t(m)] != i64(H[size_t(m)]) * words_per_row(W[size_t(m)]))
        throw std::invalid_argument("word_start does not match H * ceil(W / 32)");
      pix_start[size_t(m + 1)] = pix_start[size_t(m)] + i64(H[size_t(m)]) * W[size_t(m)];
      row_start[size_t(m + 1)] = row_start[size_t(m)] + H[size_t(m)];
      col_start[size_t(m + 1)] = col_start[size_t(m)] + W[size_t(m)];
    }
    words = word_start[size_t(n)] - word_start[0], pix = pix_start[size_t(n)], rows = row_start[size_t(n)],
    cols = col_start[size_t(n)];
  }
};

struct FrameOffsets {
  size_t words, word_start, pix_start, row_start, col_start, H, W;
};

FrameOffsets place(Layout& l, const HostFrames& h) {
  FrameOffsets o{};
  o.words = l.put<u32>(h.words), o.word_start = l.put<i64>(h.n + 1), o.pix_start = l.put<i64>(h.n + 1);
  o.row_start = l.put<i64>(h.n + 1), o.col_start = l.put<i64>(h.n + 1), o.H = l.put<i32>(h.n), o.W = l.put<i32>(h.n);
  return o;
}

Frames upload(u8* d, const FrameOffsets& o, const HostFrames& h, const k2::Bits& b, Stream st) {
  std::vector<i64> ws(h.word_start);
  for (auto& x : ws) x -= h.word_start[0];
  to_device(d + o.words, b.words + h.word_start[0], sizeof(u32) * size_t(h.words), st);
  put(d, o.word_start, ws, st), put(d, o.pix_start, h.pix_start, st), put(d, o.row_start, h.row_start, st);
  put(d, o.col_start, h.col_start, st), put(d, o.H, h.H, st), put(d, o.W, h.W, st);
  return Frames{ptr<u32>(d, o.words), ptr<i64>(d, o.word_start), ptr<i64>(d, o.pix_start), ptr<i64>(d, o.row_start),
                ptr<i64>(d, o.col_start), ptr<i32>(d, o.H), ptr<i32>(d, o.W), h.n};
}

Block<true>& arena(int which) {  // grow-only device arenas, never freed (as gpucheck's)
  static Block<true>* a = new Block<true>[4];
  return a[which];
}

}  // namespace
}  // namespace gpucheck

namespace k2 {
using namespace gpucheck;

bool cuda_build() { return gpucheck::cuda_build(); }

// ===================================================================== validity
Validity validity(const Polygons& in, i64 grid_min_edges, i64 budget_per_edge) {
  clear_error();
  const Stream st = k2_stream();
  const i64 n = in.n, nr = n ? in.shape_start[n] : 0, nv = nr ? in.ring_start[nr] : 0;
  Validity out;
  if (n == 0) return out;
  std::vector<i64> cell_start(size_t(n + 1), 0);
  std::vector<i32> g(size_t(n), 0);
  for (i64 s = 0; s < n; ++s) {
    const i64 E = in.ring_start[in.shape_start[s + 1]] - in.ring_start[in.shape_start[s]];
    i32 gs = 0;
    if (E >= grid_min_edges && E > 0) {
      gs = i32(std::ceil(std::sqrt(double(E))));
      gs = gs < 1 ? 1 : (gs > 128 ? 128 : gs);
    }
    g[size_t(s)] = gs;
    cell_start[size_t(s + 1)] = cell_start[size_t(s)] + i64(gs) * gs;
  }
  const i64 ncells = cell_start[size_t(n)];
  Layout l;
  const size_t o_xy = l.put<i32>(2 * nv), o_roles = l.put<i32>(nr), o_W = l.put<i32>(n), o_H = l.put<i32>(n),
               o_rs = l.put<i64>(nr + 1), o_ss = l.put<i64>(n + 1), o_zero = l.put<i64>(n + 1),
               o_rw = l.put<u8>(n), o_vw = l.put<u8>(n);
  size_t o_cxy[2], o_len[2], o_fewer[2], o_area2[2], o_box[2], o_counts[2], o_pairs[2], o_def[2];
  for (int k = 0; k < 2; ++k) {
    o_cxy[k] = l.put<i32>(2 * nv), o_len[k] = l.put<i32>(nr), o_fewer[k] = l.put<u8>(nr);
    o_area2[k] = l.put<i64>(nr), o_box[k] = l.put<i32>(4 * nr), o_counts[k] = l.put<i32>(nr * kNumDefects);
    o_pairs[k] = l.put<i32>(n * 2 * kPairs), o_def[k] = l.put<i64>(n * kNumDefects);
  }
  const size_t o_clean = l.put<u8>(n), o_bad = l.put<i32>(n), o_g = l.put<i32>(n), o_cs = l.put<i64>(n + 1),
               o_gbox = l.put<i32>(4 * n), o_cnt = l.put<i32>(ncells), o_off = l.put<i64>(ncells + 1),
               o_she = l.put<i64>(n), o_shb = l.put<i64>(n + 1), o_use = l.put<u8>(n);
  u8* d = arena(0).get(l.size);
  zero(d, l.size, st);
  to_device(d + o_xy, in.xy, sizeof(i32) * size_t(2 * nv), st);
  to_device(d + o_roles, in.roles, sizeof(i32) * size_t(nr), st);
  to_device(d + o_W, in.W, sizeof(i32) * size_t(n), st);
  to_device(d + o_H, in.H, sizeof(i32) * size_t(n), st);
  to_device(d + o_rs, in.ring_start, sizeof(i64) * size_t(nr + 1), st);
  to_device(d + o_ss, in.shape_start, sizeof(i64) * size_t(n + 1), st);
  put(d, o_g, g, st), put(d, o_cs, cell_start, st);
  const Shapes s{ptr<i32>(d, o_xy), ptr<i64>(d, o_rs), ptr<i64>(d, o_ss), ptr<i64>(d, o_zero), ptr<i64>(d, o_zero),
                 ptr<i32>(d, o_roles), ptr<i32>(d, o_W), ptr<i32>(d, o_H), n, in.G};
  Rings r[2];
  for (int k = 0; k < 2; ++k)
    r[k] = Rings{ptr<i32>(d, o_cxy[k]), ptr<i32>(d, o_len[k]), ptr<u8>(d, o_fewer[k]), ptr<i64>(d, o_area2[k]),
                 ptr<i32>(d, o_box[k]), ptr<i32>(d, o_counts[k]), ptr<i32>(d, o_pairs[k])};
  u8* vw = ptr<u8>(d, o_vw);
  for_each(st, n, Precision{s, ptr<u8>(d, o_rw), vw});
  for_each(st, nr, CleanRings{s, r[0]});
  for_each(st, nr, CleanRings{s, r[1]});
  Timer tm(st);
  // old: gpucheck's two steps
  tm.start(), for_each(st, nv, EdgePairs{s, r[0], vw}), out.t_old.add("EdgePairs", tm.stop());
  tm.start(), for_each(st, n, ShapeTotals{s, r[0], vw, ptr<i64>(d, o_def[0])}), out.t_old.add("ShapeTotals", tm.stop());
  // new: grid pairs (one host round trip for the entry counts), then the warp totals
  GridView gv{ptr<i32>(d, o_g), ptr<i64>(d, o_cs), ptr<i32>(d, o_gbox), ptr<i32>(d, o_cnt), ptr<i64>(d, o_off),
              ptr<i64>(d, o_she), ptr<i64>(d, o_shb), ptr<u8>(d, o_use), nullptr, nullptr, ncells};
  tm.start();
  for_each(st, n, GridBox{s, r[1], gv});
  for_each(st, nv, GridCount{s, r[1], gv});
  out.t_new.add("GridBox+Count", tm.stop());
  tm.start();
  const std::vector<i64> she = get<i64>(d, o_she, n, st);
  std::vector<i64> base(size_t(n + 1), 0);
  std::vector<u8> use(size_t(n), 0);
  for (i64 k = 0; k < n; ++k) {
    const i64 E = in.ring_start[in.shape_start[k + 1]] - in.ring_start[in.shape_start[k]];
    use[size_t(k)] = g[size_t(k)] > 0 && she[size_t(k)] <= budget_per_edge * E + 1024;
    base[size_t(k + 1)] = base[size_t(k)] + (use[size_t(k)] ? she[size_t(k)] : 0);
  }
  const i64 total = base[size_t(n)];
  put(d, o_shb, base, st), put(d, o_use, use, st);
  to_device(d + o_off + sizeof(i64) * size_t(ncells), &total, sizeof(i64), st);
  const size_t ent = sizeof(i64) * size_t(total > 0 ? total : 1);
  u8* e = arena(1).get(ent + sizeof(i32) * size_t(total > 0 ? total : 1));
  gv.entries = ptr<i64>(e, 0), gv.ering = ptr<i32>(e, ent);
  out.t_new.add("host round trip", tm.stop());
  tm.start();
  for_each(st, n, GridScan{s, gv});
  for_each(st, nv, GridFill{s, r[1], gv});
  out.t_new.add("GridScan+Fill", tm.stop());
  tm.start(), for_each(st, total, GridPairs{s, r[1], gv, vw}), out.t_new.add("GridPairs", tm.stop());
  tm.start(), for_each(st, nv, EdgePairsRest{s, r[1], vw, gv.use}), out.t_new.add("EdgePairsRest", tm.stop());
  tm.start();
  for_each_warp(st, n, ShapeTotalsWarp{s, r[1], vw, ptr<i64>(d, o_def[1]), ptr<u8>(d, o_clean), ptr<i32>(d, o_bad)});
  out.t_new.add("ShapeTotalsWarp", tm.stop());
  out.defects_old = get<i64>(d, o_def[0], n * kNumDefects, st);
  out.defects_new = get<i64>(d, o_def[1], n * kNumDefects, st);
  out.counts_old = get<i32>(d, o_counts[0], nr * kNumDefects, st);
  out.counts_new = get<i32>(d, o_counts[1], nr * kNumDefects, st);
  out.pairs_old = get<i32>(d, o_pairs[0], n * 2 * kPairs, st);
  out.pairs_new = get<i32>(d, o_pairs[1], n * 2 * kPairs, st);
  out.clean = get<u8>(d, o_clean, n, st);
  out.bad_rings = get<i32>(d, o_bad, n, st);
  out.used_grid = use;
  out.grid_entries = total;
  return out;
}

// ===================================================================== boundary, EDT
std::vector<u32> boundary(const Bits& m, std::vector<i64>* count) {
  clear_error();
  const Stream st = k2_stream();
  const HostFrames h(m);
  Layout l;
  const FrameOffsets fo = place(l, h);
  const size_t o_out = l.put<u32>(h.words), o_nb = l.put<i64>(h.n);
  u8* d = arena(0).get(l.size);
  zero(d + o_nb, sizeof(i64) * size_t(h.n), st);
  const Frames f = upload(d, fo, h, m, st);
  for_each(st, h.rows, BoundaryRows{f, ptr<u32>(d, o_out), ptr<i64>(d, o_nb)});
  if (count) *count = get<i64>(d, o_nb, h.n, st);
  return get<u32>(d, o_out, h.words, st);
}

std::vector<i32> edt_sq(const Bits& src, Times* t) {
  clear_error();
  const Stream st = k2_stream();
  const HostFrames h(src);
  Layout l;
  const FrameOffsets fo = place(l, h);
  const size_t o_g = l.put<i32>(h.pix), o_s = l.put<i32>(h.pix), o_t = l.put<i32>(h.pix), o_d = l.put<i32>(h.pix);
  u8* d = arena(0).get(l.size);
  const Frames f = upload(d, fo, h, src, st);
  Timer tm(st);
  tm.start(), for_each_warp(st, h.rows, EdtRowsH{f, ptr<i32>(d, o_g)});
  if (t) t->add("EdtRowsH", tm.stop());
  tm.start(), for_each(st, h.cols, EdtColsV{f, ptr<i32>(d, o_g), ptr<i32>(d, o_s), ptr<i32>(d, o_t), ptr<i32>(d, o_d)});
  if (t) t->add("EdtColsV", tm.stop());
  return get<i32>(d, o_d, h.pix, st);
}

// ===================================================================== boundary F, HD95
Surface surface(const Bits& m, const PairQuery& q) {
  clear_error();
  const Stream st = k2_stream();
  const HostFrames h(m);
  Surface out;
  const i64 np = q.n, nd = 2 * np;
  for (i64 p = 0; p < np; ++p) {
    if (q.a[p] < 0 || q.a[p] >= h.n || q.b[p] < 0 || q.b[p] >= h.n) throw std::invalid_argument("pair index");
    if (h.H[size_t(q.a[p])] != h.H[size_t(q.b[p])] || h.W[size_t(q.a[p])] != h.W[size_t(q.b[p])])
      throw std::invalid_argument("a pair's masks must share one frame");
  }
  if (q.tau2 < 0 || q.tau2 >= kNoSource || q.p_den <= 0 || q.p_num < 0 || q.p_num > q.p_den)
    throw std::invalid_argument("tau2 in [0, 2^31 - 1), 0 <= p_num <= p_den");
  Layout l;
  const FrameOffsets fo = place(l, h);
  const size_t o_bnd = l.put<u32>(h.words), o_nb = l.put<i64>(h.n), o_g = l.put<i32>(h.pix), o_s = l.put<i32>(h.pix),
               o_t = l.put<i32>(h.pix), o_d = l.put<i32>(h.pix), o_a = l.put<i64>(np), o_b = l.put<i64>(np),
               o_ls = l.put<i64>(nd + 1), o_rs = l.put<i64>(nd + 1), o_fill = l.put<i32>(nd),
               o_hit = l.put<i64>(nd), o_rank = l.put<i64>(nd), o_lo = l.put<i64>(nd), o_hi = l.put<i64>(nd);
  u8* d = arena(0).get(l.size);
  zero(d + o_nb, sizeof(i64) * size_t(h.n), st);
  zero(d + o_fill, sizeof(i32) * size_t(nd), st);
  zero(d + o_hit, sizeof(i64) * size_t(nd), st);
  const Frames f = upload(d, fo, h, m, st);
  Timer tm(st);
  tm.start(), for_each(st, h.rows, BoundaryRows{f, ptr<u32>(d, o_bnd), ptr<i64>(d, o_nb)});
  out.t.add("BoundaryRows", tm.stop());
  const Frames fb{ptr<u32>(d, o_bnd), f.word_start, f.pix_start, f.row_start, f.col_start, f.H, f.W, f.n};
  tm.start(), for_each_warp(st, h.rows, EdtRowsH{fb, ptr<i32>(d, o_g)}), out.t.add("EdtRowsH", tm.stop());
  tm.start(), for_each(st, h.cols, EdtColsV{fb, ptr<i32>(d, o_g), ptr<i32>(d, o_s), ptr<i32>(d, o_t), ptr<i32>(d, o_d)});
  out.t.add("EdtColsV", tm.stop());
  tm.start();
  const std::vector<i64> nb = get<i64>(d, o_nb, h.n, st);  // round trip: the list sizes
  std::vector<i64> ls(size_t(nd + 1), 0), rs(size_t(nd + 1), 0), A(q.a, q.a + np), B(q.b, q.b + np);
  for (i64 k = 0; k < nd; ++k) {
    const i64 src = (k & 1) ? q.b[k >> 1] : q.a[k >> 1];
    ls[size_t(k + 1)] = ls[size_t(k)] + nb[size_t(src)];
    rs[size_t(k + 1)] = rs[size_t(k)] + h.H[size_t(src)];
  }
  put(d, o_a, A, st), put(d, o_b, B, st), put(d, o_ls, ls, st), put(d, o_rs, rs, st);
  u8* lst = arena(1).get(sizeof(i32) * size_t(ls[size_t(nd)] > 0 ? ls[size_t(nd)] : 1));
  out.t.add("host round trip", tm.stop());
  const Dirs dirs{ptr<i64>(d, o_a), ptr<i64>(d, o_b), ptr<i64>(d, o_ls), ptr<i64>(d, o_rs), nd};
  tm.start();
  for_each(st, rs[size_t(nd)], Directed{f, dirs, ptr<u32>(d, o_bnd), ptr<i32>(d, o_d), q.tau2, ptr<i32>(d, o_fill),
                                        ptr<i64>(d, o_hit), ptr<i32>(lst, 0)});
  out.t.add("Directed", tm.stop());
  tm.start();
  for_each_warp(st, nd, Select{dirs, ptr<i32>(lst, 0), q.p_num, q.p_den, ptr<i64>(d, o_rank), ptr<i64>(d, o_lo),
                               ptr<i64>(d, o_hi)});
  out.t.add("Select", tm.stop());
  const std::vector<i64> hit = get<i64>(d, o_hit, nd, st), rank = get<i64>(d, o_rank, nd, st),
                         lo = get<i64>(d, o_lo, nd, st), hi = get<i64>(d, o_hi, nd, st);
  for (i64 p = 0; p < np; ++p) {
    const size_t a = size_t(2 * p), b = a + 1;
    out.n_a.push_back(ls[a + 1] - ls[a]), out.n_b.push_back(ls[b + 1] - ls[b]);
    out.hit_a.push_back(hit[a]), out.hit_b.push_back(hit[b]);
    out.rank_a.push_back(rank[a]), out.rank_b.push_back(rank[b]);
    out.lo_a.push_back(lo[a]), out.hi_a.push_back(hi[a]), out.lo_b.push_back(lo[b]), out.hi_b.push_back(hi[b]);
  }
  return out;
}

double boundary_f(i64 n_a, i64 n_b, i64 hit_a, i64 hit_b) {
  if (n_a == 0 && n_b == 0) return 1.0;
  if (n_a == 0 || n_b == 0) return 0.0;
  const double p = double(hit_a) / double(n_a), r = double(hit_b) / double(n_b);
  return p + r > 0 ? 2.0 * p * r / (p + r) : 0.0;
}

static double directed_percentile(i64 n, i64 rank, i64 lo, i64 hi, i64 p_num, i64 p_den) {
  if (n == 0 || lo == kNoSource || hi == kNoSource) return std::numeric_limits<double>::infinity();
  const double dl = std::sqrt(double(lo)), dh = std::sqrt(double(hi));
  const double frac = double(p_num * (n - 1) - rank * p_den) / double(p_den);
  return dl + (dh - dl) * frac;
}

double hd_percentile(i64 n_a, i64 rank_a, i64 lo_a, i64 hi_a, i64 n_b, i64 rank_b, i64 lo_b, i64 hi_b, i64 p_num,
                     i64 p_den) {
  const double x = directed_percentile(n_a, rank_a, lo_a, hi_a, p_num, p_den);
  const double y = directed_percentile(n_b, rank_b, lo_b, hi_b, p_num, p_den);
  return x > y ? x : y;
}

// ===================================================================== sets
SetStats sets(const Bits& m, const SetQuery& q) {
  clear_error();
  const Stream st = k2_stream();
  const HostFrames h(m);
  SetStats out;
  const i64 n = q.n;
  std::vector<i64> pair_start(size_t(n + 1), 0), state_start(size_t(n + 1), 0);
  out.inter_start.assign(size_t(n + 1), 0), out.choice_start.assign(size_t(n + 1), 0);
  for (i64 s = 0; s < n; ++s) {
    const i64 m0 = q.mask_start[s], M = q.mask_start[s + 1] - m0, K = q.K[s];
    if (m0 < 0 || q.mask_start[s + 1] > h.n || K < 1 || M - K < 1) throw std::invalid_argument("set: K >= 1, k >= 1");
    for (i64 j = m0; j < m0 + M; ++j)
      if (h.H[size_t(j)] != h.H[size_t(m0)] || h.W[size_t(j)] != h.W[size_t(m0)])
        throw std::invalid_argument("a set's masks must share one frame");
    pair_start[size_t(s + 1)] = pair_start[size_t(s)] + M * (M + 1) / 2;
    out.inter_start[size_t(s + 1)] = out.inter_start[size_t(s)] + M * M;
    out.choice_start[size_t(s + 1)] = out.choice_start[size_t(s)] + K;
    state_start[size_t(s + 1)] = state_start[size_t(s)] + hm_states(K, M - K);
  }
  std::vector<i64> ms(q.mask_start, q.mask_start + n + 1), ws(h.word_start);
  for (auto& x : ws) x -= h.word_start[0];
  const std::vector<i32> K(q.K, q.K + n);
  Layout l;
  const size_t o_words = l.put<u32>(h.words), o_ws = l.put<i64>(h.n + 1), o_H = l.put<i32>(h.n),
               o_W = l.put<i32>(h.n), o_ms = l.put<i64>(n + 1), o_K = l.put<i32>(n), o_ps = l.put<i64>(n + 1),
               o_is = l.put<i64>(n + 1), o_ss = l.put<i64>(n + 1), o_cs = l.put<i64>(n + 1),
               o_inter = l.put<i64>(out.inter_start[size_t(n)]), o_best = l.put<u64>(state_start[size_t(n)]),
               o_choice = l.put<i32>(out.choice_start[size_t(n)]), o_hq = l.put<u64>(n), o_hm = l.put<double>(n),
               o_ged = l.put<double>(n);
  u8* d = arena(0).get(l.size);
  to_device(d + o_words, m.words + h.word_start[0], sizeof(u32) * size_t(h.words), st);
  put(d, o_ws, ws, st), put(d, o_H, h.H, st), put(d, o_W, h.W, st), put(d, o_ms, ms, st), put(d, o_K, K, st);
  put(d, o_ps, pair_start, st), put(d, o_is, out.inter_start, st), put(d, o_ss, state_start, st);
  put(d, o_cs, out.choice_start, st);
  const SetView v{ptr<u32>(d, o_words), ptr<i64>(d, o_ws), ptr<i32>(d, o_H), ptr<i32>(d, o_W), ptr<i64>(d, o_ms),
                  ptr<i32>(d, o_K), ptr<i64>(d, o_ps), ptr<i64>(d, o_is), n};
  Timer tm(st);
  tm.start(), for_each_warp(st, pair_start[size_t(n)], SetInter{v, ptr<i64>(d, o_inter)});
  out.t.add("SetInter", tm.stop());
  tm.start();
  for_each(st, n, SetTerms{v, ptr<i64>(d, o_inter), ptr<i64>(d, o_ss), ptr<i64>(d, o_cs), ptr<u64>(d, o_best),
                           ptr<i32>(d, o_choice), ptr<u64>(d, o_hq), ptr<double>(d, o_hm), ptr<double>(d, o_ged)});
  out.t.add("SetTerms", tm.stop());
  out.inter = get<i64>(d, o_inter, out.inter_start[size_t(n)], st);
  out.choice = get<i32>(d, o_choice, out.choice_start[size_t(n)], st);
  out.hm_q = get<u64>(d, o_hq, n, st);
  out.hm_iou = get<double>(d, o_hm, n, st);
  out.ged = get<double>(d, o_ged, n, st);
  return out;
}

}  // namespace k2
