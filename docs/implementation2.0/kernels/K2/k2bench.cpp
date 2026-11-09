// k2bench.cpp -- micro-benchmark of the K2 kernels. Per workload: device ms per stage (CUDA events on the GPU build,
// steady clock on the CPU build; best of `reps`) and the wall time of the whole API call (copies included).
// Usage: k2bench [scale] [reps]   scale 1 = the GPU run (512 x 512 frames), scale 8 = 8x smaller batches (CPU build)
#include <chrono>
#include <cstdio>
#include <cstdlib>

#include "k2.h"
#include "k2gen.h"

using namespace k2gen;

static double now_ms() {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static void report(const char* what, const k2::Times& best, double wall, const char* unit, double items) {
  std::printf("%-34s wall %9.3f ms  device %9.3f ms  (%.1f %s/s device)\n", what, wall, best.total(),
              items / (best.total() / 1000.0), unit);
  for (int i = 0; i < best.n; ++i) std::printf("    %-22s %9.3f ms\n", best.name[i], best.ms[i]);
}

template <class F>
static void bench(const char* what, int reps, const char* unit, double items, F f) {
  k2::Times best;
  double wall = 1e300;
  for (int r = 0; r < reps; ++r) {
    const double t0 = now_ms();
    const k2::Times t = f();
    const double w = now_ms() - t0;
    if (w < wall) wall = w;
    if (r == 0 || t.total() < best.total()) best = t;
  }
  report(what, best, wall, unit, items);
}

int main(int argc, char** argv) {
  const int scale = argc > 1 ? std::atoi(argv[1]) : 1, reps = argc > 2 ? std::atoi(argv[2]) : 3;
  std::printf("k2bench: %s build, scale 1/%d, best of %d\n", k2::cuda_build() ? "CUDA" : "CPU", scale, reps);
  Rng r(424242);
  const i32 F = 512;  // LIDC frame

  // Items: K = 8 samples (jittered rater copies) + k = 4 raters on one 512 x 512 frame.
  const int items = 32 / scale > 0 ? 32 / scale : 1, K = 8, k = 4;
  std::vector<Mask> ms;
  std::vector<i64> start{0}, A, B;
  std::vector<i32> Ks;
  for (int it = 0; it < items; ++it) {
    Mask base = blank(F, F);
    const i64 cx = r.range(150, 350), cy = r.range(150, 350), rad = r.range(10, 90);
    for (i32 y = 0; y < F; ++y)
      for (i32 x = 0; x < F; ++x)
        if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= rad * rad) base.px[size_t(y) * F + x] = 1;
    const i64 first = i64(ms.size());
    for (int i = 0; i < K + k; ++i) ms.push_back(jitter(r, base, 30));
    for (int i = 0; i < K; ++i)
      for (int j = 0; j < k; ++j) A.push_back(first + i), B.push_back(first + K + j);
    start.push_back(i64(ms.size()));
    Ks.push_back(K);
  }
  const Packed P(ms);
  const double npx = double(ms.size()) * F * F;

  bench("edt_sq (all masks)", reps, "Mpx", npx / 1e6, [&] {
    k2::Times t;
    (void)k2::edt_sq(P.bits(), &t);
    return t;
  });
  bench("surface: boundary F + HD95", reps, "pairs", double(A.size()), [&] {
    k2::PairQuery q;
    q.n = i64(A.size()), q.a = A.data(), q.b = B.data(), q.tau2 = 4;
    return k2::surface(P.bits(), q).t;
  });
  bench("sets: GED + HM-IoU (K=8, k=4)", reps, "sets", double(items), [&] {
    k2::SetQuery q;
    q.n = i64(Ks.size()), q.mask_start = start.data(), q.K = Ks.data();
    return k2::sets(P.bits(), q).t;
  });

  // Validity, typical: one star ring of 64-200 vertices per shape (model outlines on the 0-1000 grid).
  {
    Polys T;
    const int n = 4096 / scale;
    for (int i = 0; i < n; ++i) {
      T.ring(star_ring(r, int(r.range(64, 200)), 500, 500, 50, 400, false), -1);
      T.end_shape(512, 512);
    }
    const gpucheck::Polygons in = T.view();
    for (const i64 min_edges : {i64(1) << 30, i64(0)}) {
      k2::Validity best;
      double wall = 1e300;
      for (int rep = 0; rep < reps; ++rep) {
        const double t0 = now_ms();
        k2::Validity v = k2::validity(in, min_edges, 64);
        wall = std::min(wall, now_ms() - t0);
        if (rep == 0 || v.t_new.total() < best.t_new.total()) best = std::move(v);
      }
      std::printf("validity typical (%d shapes), grid %s, wall (both paths) %.3f ms\n", n,
                  min_edges ? "off" : "on", wall);
      report("  old: EdgePairs + ShapeTotals", best.t_old, 0, "shapes", n);
      report("  new: grid / rest + warp totals", best.t_new, 0, "shapes", n);
    }
  }
  // Validity, adversarial: strips (touching holes: the witness straggler) and 1000-vertex stars with holes.
  {
    Polys T;
    const int n = 32 / scale > 0 ? 32 / scale : 1;
    for (int i = 0; i < n; ++i) strips_shape(T, 1000, 9, 1000);
    for (int i = 0; i < n; ++i) {
      T.ring(star_ring(r, 1000, 500, 500, 300, 499, false), -1);
      for (int h = 0; h < 9; ++h)
        T.ring(star_ring(r, 1000, 300 + 50 * (h % 3), 300 + 50 * (h / 3), 1, 20, true), 0);
      T.end_shape(512, 512);
    }
    const gpucheck::Polygons in = T.view();
    k2::Validity best;
    double wall = 1e300;
    for (int rep = 0; rep < reps; ++rep) {
      const double t0 = now_ms();
      k2::Validity v = k2::validity(in, 512, 64);
      wall = std::min(wall, now_ms() - t0);
      if (rep == 0 || v.t_new.total() < best.t_new.total()) best = std::move(v);
    }
    i64 used = 0;
    for (auto u : best.used_grid) used += u;
    std::printf("validity adversarial (%d shapes x 10 rings x 1000 vertices, grid on %lld), wall %.3f ms, same=%d\n",
                2 * n, (long long)used, wall, int(best.defects_new == best.defects_old));
    report("  old: EdgePairs + ShapeTotals", best.t_old, 0, "shapes", 2 * n);
    report("  new: grid / rest + warp totals", best.t_new, 0, "shapes", 2 * n);
  }
  return 0;
}
