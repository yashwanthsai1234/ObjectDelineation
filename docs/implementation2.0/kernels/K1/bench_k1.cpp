// bench_k1.cpp -- K1 versus gpucheck::check_masks (the dense slot-grid polygonizer it replaces), same build type:
// both CPU (clang++/g++) or both CUDA (nvcc). Every timed result is also checked for exact equality.
//   usage: bench_k1 [reps=10] [real_prefix]       real_prefix: raw masks written by cross_numpy.py (<prefix>.bits ...)
// Output: one line per workload: median ms per call, masks/s, K1 stage split (h2d / ccl / rings / trace / d2h, ms).
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "gpucheck.h"
#include "k1poly.h"
#include "masks_gen.h"

using namespace k1gen;

template <class T>
static std::vector<T> slurp(const std::string& path) {
  std::ifstream f(path, std::ios::binary | std::ios::ate);
  if (!f) return {};
  const auto n = size_t(f.tellg());
  std::vector<T> v(n / sizeof(T));
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), std::streamsize(n));
  return v;
}

static double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return v.empty() ? 0 : v[v.size() / 2];
}

static double time_ms(const std::function<void()>& f) {
  const auto t0 = std::chrono::steady_clock::now();
  f();
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

static void run(const char* name, const Batch& b, int reps, int seed_rows) {
  const auto in = b.view();
  const gpucheck::Masks gin{in.n, in.bits, in.start, in.H, in.W};
  k1poly::Options opt;
  opt.seed_rows = seed_rows;
  k1poly::Result r;
  gpucheck::MaskResult g;
  for (int w = 0; w < 2; ++w) r = k1poly::polygonize_masks(in, opt), g = gpucheck::check_masks(gin);  // warm-up
  const bool eq = g.quads == r.quads && g.ring_start == r.ring_start && g.pts_start == r.pts_start &&
                  g.ring_area == r.ring_area && g.pts == r.pts;
  std::vector<double> tk, tg;
  k1poly::Timing sum;
  for (int k = 0; k < reps; ++k) {
    k1poly::Timing t;
    tk.push_back(time_ms([&] { r = k1poly::polygonize_masks(in, opt, &t); }));
    sum.h2d += t.h2d, sum.ccl += t.ccl, sum.rings += t.rings, sum.trace += t.trace, sum.d2h += t.d2h;
    tg.push_back(time_ms([&] { g = gpucheck::check_masks(gin); }));
  }
  const double mk = median(tk), mg = median(tg);
  std::printf("%-26s n=%-5lld seeds=%-2d rings=%-8zu corners=%-9zu | K1 %9.3f ms (%9.0f masks/s) | check_masks %9.3f ms "
              "(%9.0f masks/s) | x%.1f | %s | K1 split h2d %.3f ccl %.3f rings %.3f trace %.3f d2h %.3f\n",
              name, (long long)in.n, seed_rows, r.ring_area.size(), r.pts.size() / 2, mk, in.n / (mk / 1e3), mg,
              in.n / (mg / 1e3), mg / std::max(mk, 1e-9), eq ? "IDENTICAL" : "DIFFERENT", sum.h2d / reps,
              sum.ccl / reps, sum.rings / reps, sum.trace / reps, sum.d2h / reps);
  std::fflush(stdout);
}

int main(int argc, char** argv) {
  const int reps = argc > 1 ? std::atoi(argv[1]) : 10;
  std::printf("bench_k1: %s build (k1poly %s, gpucheck %s), reps %d\n", k1poly::cuda_build() ? "CUDA" : "CPU",
              k1poly::cuda_build() ? "cuda" : "cpu", gpucheck::cuda_build() ? "cuda" : "cpu", reps);
  std::printf("device: %s\n", k1poly::device_info().c_str());
  std::mt19937_64 rng(5);
  auto many = [&](int n, const std::function<Mask()>& make) {
    Batch b;
    for (int k = 0; k < n; ++k) b.add(make());
    return b;
  };
  // 1. lesion-like: 512 x 512 frames with a few blobs (the RL / eval frame size), batch 256
  const Batch lesion = many(256, [&] { return blobs(rng, 512, 512, 2, 18.0); });
  run("lesion512 b256", lesion, reps, 16);
  run("lesion512 b256", lesion, reps, 0);
  // 2. many rings: noise p = 0.5, 512 x 512, batch 16 (~ 25 k rings per mask)
  const Batch noisy = many(16, [&] { return noise(rng, 512, 512, 0.5); });
  run("noise512 b16", noisy, reps, 16);
  // 3. one very long ring (comb, ~260 k edges): seed rows split its trace; 0 = one thread per ring
  const Batch combs = many(4, [&] { return comb(512, 512); });
  run("comb512 b4", combs, reps, 16);
  run("comb512 b4", combs, reps, 0);
  // 4. deep union-find chains and nested rings
  const Batch spirals = many(4, [&] { return spiral(512, 512); });
  run("spiral512 b4", spirals, reps, 16);
  const Batch conc = many(4, [&] { return concentric(512, 512); });
  run("concentric512 b4", conc, reps, 16);
  const Batch chk = many(4, [&] { return checker(512, 512); });
  run("checker512 b4", chk, reps, 16);
  // 5. real masks (LIDC crops of data/checker_eval/masks.npz, exported by cross_numpy.py), if given
  if (argc > 2) {
    const std::string pre = argv[2];
    auto bits = slurp<u8>(pre + ".bits");
    auto hw = slurp<i32>(pre + ".hw");
    auto start = slurp<i64>(pre + ".start");
    if (!hw.empty()) {
      Batch real;
      real.bits = bits, real.start = start;
      for (size_t k = 0; k < hw.size() / 2; ++k) real.H.push_back(hw[2 * k]), real.W.push_back(hw[2 * k + 1]);
      run("real (all)", real, reps, 16);
    }
  }
  return 0;
}
