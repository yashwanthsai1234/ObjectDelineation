// k2test.cpp -- exact-equality tests of the K2 kernels against independent references, on random and adversarial
// inputs, plus one FNV-1a digest of every output: the CPU build and the GPU build must print the same digest for the
// same seed. Usage: k2test [seed] [scale]   (scale 1 = quick, 4 = the GPU run's size)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "k2.h"
#include "k2gen.h"
#include "k2ref.h"

using namespace k2gen;

static int g_fail = 0;
static long g_checks = 0;

template <class T>
static bool same(const char* what, const std::vector<T>& a, const std::vector<T>& b, i64 context = -1) {
  ++g_checks;
  if (a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), sizeof(T) * a.size()) == 0)) return true;
  size_t i = 0;
  while (i < a.size() && i < b.size() && std::memcmp(&a[i], &b[i], sizeof(T)) == 0) ++i;
  std::printf("FAIL %s (context %lld): sizes %zu vs %zu, first difference at %zu\n", what, (long long)context, a.size(),
              b.size(), i);
  ++g_fail;
  return false;
}

static bool same_i(const char* what, i64 a, i64 b, i64 context) {
  ++g_checks;
  if (a == b) return true;
  std::printf("FAIL %s (context %lld): %lld vs %lld\n", what, (long long)context, (long long)a, (long long)b);
  ++g_fail;
  return false;
}

static bool same_bits(const char* what, double a, double b, i64 context) {  // bit-exact float64
  ++g_checks;
  if (std::memcmp(&a, &b, sizeof a) == 0) return true;
  std::printf("FAIL %s (context %lld): %.17g vs %.17g\n", what, (long long)context, a, b);
  ++g_fail;
  return false;
}

// ---------------------------------------------------------------- EDT and boundary
static void test_edt(Rng& r, Digest& dg, int scale) {
  std::vector<Mask> ms;
  const i32 dims[][2] = {{1, 1}, {1, 7}, {7, 1}, {2, 33}, {33, 2}, {31, 32}, {32, 31}, {17, 63}, {40, 64}, {33, 65},
                         {48, 48}, {5, 100}, {100, 5}};
  for (const auto& d : dims)
    for (int t = 0; t < 4 * scale; ++t) ms.push_back(random_mask(r, d[0], d[1]));
  for (int t = 0; t < 2 * scale; ++t) {  // larger, sparse sources: long distances
    Mask m = blank(i32(r.range(60, 130)), i32(r.range(60, 140)));
    for (int k = int(r.range(0, 4)); k > 0; --k) m.px[size_t(r.range(0, m.H - 1)) * m.W + size_t(r.range(0, m.W - 1))] = 1;
    ms.push_back(m);
  }
  const Packed P(ms);
  std::vector<i64> nb;
  const std::vector<u32> bw = k2::boundary(P.bits(), &nb);
  const std::vector<i32> d2 = k2::edt_sq(P.bits());
  dg.vec(bw), dg.vec(nb), dg.vec(d2);
  i64 at = 0;
  for (size_t m = 0; m < ms.size(); ++m) {
    const std::vector<i32> ref = k2ref::edt_brute(ms[m]);
    const std::vector<i32> got(d2.begin() + at, d2.begin() + at + i64(ref.size()));
    same("edt_sq vs brute force", got, ref, i64(m));
    at += i64(ref.size());
    const Mask bref = k2ref::boundary(ms[m]);
    std::vector<u32> wref, w;
    std::vector<i64> ws;
    std::vector<i32> H, W;
    k2ref::pack({bref}, wref, ws, H, W);
    w.assign(bw.begin() + P.word_start[m], bw.begin() + P.word_start[m + 1]);
    same("boundary words", w, wref, i64(m));
    i64 c = 0;
    for (auto v : bref.px) c += v;
    same_i("boundary count", nb[m], c, i64(m));
  }
  std::printf("edt+boundary: %zu masks\n", ms.size());
}

// ---------------------------------------------------------------- boundary F, HD95
static void test_surface(Rng& r, Digest& dg, int scale) {
  std::vector<Mask> ms;
  std::vector<i64> A, B;
  const i32 dims[][2] = {{1, 1}, {3, 40}, {32, 32}, {33, 47}, {48, 64}, {64, 33}};
  for (const auto& d : dims)
    for (int grp = 0; grp < 3 * scale; ++grp) {
      const i64 first = i64(ms.size());
      const Mask base = random_mask(r, d[0], d[1]);
      ms.push_back(base);
      for (int t = 0; t < 3; ++t) ms.push_back(r.chance(70) ? jitter(r, base, int(r.range(5, 60))) : random_mask(r, d[0], d[1]));
      for (i64 a = first; a < i64(ms.size()); ++a)
        for (i64 b = first; b < i64(ms.size()); ++b) A.push_back(a), B.push_back(b);
    }
  const Packed P(ms);
  const i64 tau2s[] = {0, 1, 2, 4, 8};
  const i64 pcts[][2] = {{95, 100}, {1, 2}, {100, 100}, {0, 1}};
  for (int v = 0; v < 4; ++v) {
    k2::PairQuery q;
    q.n = i64(A.size()), q.a = A.data(), q.b = B.data(), q.tau2 = tau2s[v % 5], q.p_num = pcts[v][0], q.p_den = pcts[v][1];
    const k2::Surface s = k2::surface(P.bits(), q);
    dg.vec(s.n_a), dg.vec(s.n_b), dg.vec(s.hit_a), dg.vec(s.hit_b), dg.vec(s.rank_a), dg.vec(s.rank_b);
    dg.vec(s.lo_a), dg.vec(s.hi_a), dg.vec(s.lo_b), dg.vec(s.hi_b);
    for (i64 p = 0; p < q.n; ++p) {
      const k2ref::Directed x = k2ref::directed(ms[size_t(A[p])], ms[size_t(B[p])], q.tau2, q.p_num, q.p_den);
      const k2ref::Directed y = k2ref::directed(ms[size_t(B[p])], ms[size_t(A[p])], q.tau2, q.p_num, q.p_den);
      const bool ok = same_i("n_a", s.n_a[p], x.n, p) && same_i("hit_a", s.hit_a[p], x.hit, p) &&
                      same_i("rank_a", s.rank_a[p], x.rank, p) && same_i("lo_a", s.lo_a[p], x.lo, p) &&
                      same_i("hi_a", s.hi_a[p], x.hi, p) && same_i("n_b", s.n_b[p], y.n, p) &&
                      same_i("hit_b", s.hit_b[p], y.hit, p) && same_i("rank_b", s.rank_b[p], y.rank, p) &&
                      same_i("lo_b", s.lo_b[p], y.lo, p) && same_i("hi_b", s.hi_b[p], y.hi, p);
      if (!ok) break;
    }
  }
  std::printf("surface: %zu pairs x 4 settings\n", A.size());
}

// ---------------------------------------------------------------- GED, HM-IoU
static void test_sets(Rng& r, Digest& dg, int scale) {
  std::vector<Mask> ms;
  std::vector<i64> start{0};
  std::vector<i32> K;
  for (int t = 0; t < 40 * scale; ++t) {
    const i32 H = i32(r.range(1, 48)), W = i32(r.range(1, 70));
    const int k = int(r.range(1, 4)), Ks = int(r.range(1, 8));
    const Mask base = random_mask(r, H, W);
    std::vector<Mask> raters;
    for (int j = 0; j < k; ++j) raters.push_back(r.chance(20) ? base : jitter(r, base, int(r.range(5, 50))));
    for (int i = 0; i < Ks; ++i) {  // samples: copies of raters (exact ties), jitters, empties
      const int c = int(r.range(0, 9));
      ms.push_back(c < 3 ? raters[size_t(r.range(0, k - 1))] : c < 8 ? jitter(r, base, int(r.range(5, 60))) : blank(H, W));
    }
    for (const Mask& m : raters) ms.push_back(m);
    start.push_back(i64(ms.size()));
    K.push_back(Ks);
  }
  {  // too many quota states: K = 1, k = 17 -> 2^17 > 2^16: the device declines (choice -1)
    for (int j = 0; j < 18; ++j) ms.push_back(random_mask(r, 8, 8));
    start.push_back(i64(ms.size()));
    K.push_back(1);
  }
  const Packed P(ms);
  k2::SetQuery q;
  q.n = i64(K.size()), q.mask_start = start.data(), q.K = K.data();
  const k2::SetStats s = k2::sets(P.bits(), q);
  dg.vec(s.inter), dg.vec(s.choice), dg.vec(s.hm_q), dg.vec(s.hm_iou), dg.vec(s.ged);
  for (i64 t = 0; t < q.n; ++t) {
    const std::vector<Mask> set(ms.begin() + start[size_t(t)], ms.begin() + start[size_t(t + 1)]);
    const i64 M = i64(set.size());
    const std::vector<i64> inter(s.inter.begin() + s.inter_start[size_t(t)], s.inter.begin() + s.inter_start[size_t(t + 1)]);
    const std::vector<i32> ch(s.choice.begin() + s.choice_start[size_t(t)], s.choice.begin() + s.choice_start[size_t(t + 1)]);
    if (M - K[size_t(t)] > 16) {
      same_i("declined set: choice", ch[0], -1, t);
      continue;
    }
    const k2ref::Set ref = k2ref::set_stats(set, K[size_t(t)]);
    same("inter matrix", inter, ref.inter, t);
    same("hm choice (lexicographic optimum)", ch, ref.choice, t);
    same_i("hm_q", i64(s.hm_q[size_t(t)]), i64(ref.hm_q), t);
    same_bits("hm_iou", s.hm_iou[size_t(t)], ref.hm_iou, t);
    same_bits("ged", s.ged[size_t(t)], ref.ged, t);
  }
  std::printf("sets: %lld sets\n", (long long)q.n);
}

// ---------------------------------------------------------------- validity
static void test_validity(Rng& r, Digest& dg, int scale) {
  Polys P;
  const int n = 150 * scale;
  for (int i = 0; i < n; ++i) random_shape(r, P, 1000, i % 10 == 0 ? 1000 : 120);
  strips_shape(P, 1000, 9, 400);
  strips_shape(P, 1000, 3, 1000);
  const gpucheck::Polygons in = P.view();
  const gpucheck::PolygonResult api = gpucheck::check_polygons(in);  // the existing path, end to end
  const struct { i64 min_edges, budget; } runs[] = {{0, 1 << 20}, {64, 4}, {1 << 30, 1}};  // all grid, mixed, none
  for (const auto& run : runs) {
    const k2::Validity v = k2::validity(in, run.min_edges, run.budget);
    same("defects: old vs check_polygons", v.defects_old, api.defects);
    same("defects: new vs old", v.defects_new, v.defects_old);
    same("counts: new vs old", v.counts_new, v.counts_old);
    same("pairs: new vs old", v.pairs_new, v.pairs_old);
    std::vector<unsigned char> clean_ref(size_t(in.n));
    std::vector<i32> bad_ref(size_t(in.n), 0);
    for (i64 s = 0; s < in.n; ++s) {
      bool z = true;
      for (int d = 0; d < gpucheck::kNumDefects; ++d) z = z && v.defects_old[size_t(s * gpucheck::kNumDefects + d)] == 0;
      clean_ref[size_t(s)] = z;
      for (i64 g = in.shape_start[s]; g < in.shape_start[s + 1]; ++g) {
        bool nz = false;
        for (int d = 0; d < gpucheck::kNumDefects; ++d) nz = nz || v.counts_old[size_t(g * gpucheck::kNumDefects + d)];
        bad_ref[size_t(s)] += nz;
      }
    }
    same("clean flag", v.clean, clean_ref);
    same("bad rings (ballot + popc)", v.bad_rings, bad_ref);
    if (run.min_edges == 0) {  // coverage: shapes with each defect class
      std::printf("defect coverage (shapes per class):");
      for (int d = 0; d < gpucheck::kNumDefects; ++d) {
        i64 c = 0;
        for (i64 s = 0; s < in.n; ++s) c += v.defects_old[size_t(s * gpucheck::kNumDefects + d)] != 0;
        std::printf(" %s=%lld", gpucheck::kDefectNames[d], (long long)c);
      }
      std::printf("\n");
    }
    i64 used = 0;
    for (auto u : v.used_grid) used += u;
    std::printf("validity: %lld shapes, grid on %lld, %lld entries\n", (long long)in.n, (long long)used,
                (long long)v.grid_entries);
    dg.vec(v.defects_new), dg.vec(v.counts_new), dg.vec(v.pairs_new), dg.vec(v.clean), dg.vec(v.bad_rings);
  }
}

int main(int argc, char** argv) {
  const k2::u64 seed = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 20261108ull;
  const int scale = argc > 2 ? std::atoi(argv[2]) : 1;
  std::printf("k2test: %s build, seed %llu, scale %d\n", k2::cuda_build() ? "CUDA" : "CPU", (unsigned long long)seed,
              scale);
  Rng r(seed);
  Digest dg;
  test_edt(r, dg, scale);
  test_surface(r, dg, scale);
  test_sets(r, dg, scale);
  test_validity(r, dg, scale);
  std::printf("checks %ld, failures %d\nDIGEST %016llx\n%s\n", g_checks, g_fail, (unsigned long long)dg.h,
              g_fail ? "FAILED" : "ALL PASS");
  return g_fail ? 1 : 0;
}
