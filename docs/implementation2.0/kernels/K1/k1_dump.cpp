// k1_dump.cpp -- run K1 on a packed mask batch read from raw files and write the result as raw arrays, for the
// cross-check against the numpy polygonizer (cross_numpy.py). With -DWITH_GPUCHECK it also compares against
// gpucheck::check_masks in-process.
//   usage: k1_dump <prefix> <out.bin> [seed_rows=16]
//   input : <prefix>.bits (uint8), <prefix>.hw (int32 n x 2: H, W), <prefix>.start (int64 n + 1)
//   output: int64 n, R, P; quads int64 (n, 5); ring_start int64 (n + 1); pts_start int64 (R + 1);
//           ring_area int64 (R); pts int32 (P, 2)
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "k1poly.h"
#ifdef WITH_GPUCHECK
#include "gpucheck.h"
#endif

using namespace k1poly;

template <class T>
static std::vector<T> slurp(const std::string& path) {
  std::ifstream f(path, std::ios::binary | std::ios::ate);
  if (!f) throw std::runtime_error("cannot open " + path);
  const auto n = size_t(f.tellg());
  if (n % sizeof(T)) throw std::runtime_error("bad size " + path);
  std::vector<T> v(n / sizeof(T));
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), std::streamsize(n));
  return v;
}

template <class T>
static void put(std::ofstream& f, const std::vector<T>& v) {
  f.write(reinterpret_cast<const char*>(v.data()), std::streamsize(v.size() * sizeof(T)));
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: k1_dump <prefix> <out.bin> [seed_rows]\n");
    return 2;
  }
  const std::string pre = argv[1];
  const auto bits = slurp<u8>(pre + ".bits");
  const auto hw = slurp<i32>(pre + ".hw");
  const auto start = slurp<i64>(pre + ".start");
  const i64 n = i64(hw.size() / 2);
  if (i64(start.size()) != n + 1 || start[n] != i64(bits.size())) throw std::runtime_error("inconsistent input");
  std::vector<i32> H(n), W(n);
  for (i64 k = 0; k < n; ++k) H[k] = hw[2 * k], W[k] = hw[2 * k + 1];
  Masks in{n, bits.data(), start.data(), H.data(), W.data()};
  Options opt;
  if (argc > 3) opt.seed_rows = std::atoi(argv[3]);
  const auto t0 = std::chrono::steady_clock::now();
  const Result r = polygonize_masks(in, opt);
  const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  const i64 R = i64(r.ring_area.size()), P = i64(r.pts.size() / 2);
  std::fprintf(stderr, "k1_dump: %lld masks, %lld rings, %lld corners, %.1f ms (%s build)\n", (long long)n,
               (long long)R, (long long)P, ms, cuda_build() ? "CUDA" : "CPU");
#ifdef WITH_GPUCHECK
  const gpucheck::MaskResult g = gpucheck::check_masks({n, bits.data(), start.data(), H.data(), W.data()});
  const bool eq = g.quads == r.quads && g.ring_start == r.ring_start && g.pts_start == r.pts_start &&
                  g.ring_area == r.ring_area && g.pts == r.pts;
  std::fprintf(stderr, "k1_dump: gpucheck::check_masks %s\n", eq ? "IDENTICAL" : "DIFFERS");
  if (!eq) return 1;
#endif
  std::ofstream f(argv[2], std::ios::binary);
  const std::vector<i64> head{n, R, P};
  put(f, head), put(f, r.quads), put(f, r.ring_start), put(f, r.pts_start), put(f, r.ring_area), put(f, r.pts);
  return 0;
}
