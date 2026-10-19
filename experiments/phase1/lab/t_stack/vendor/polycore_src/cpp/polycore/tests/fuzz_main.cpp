// Native fuzz / differential driver for polycore.hpp. Build with sanitizers:
//   clang++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all
//           -I cpp cpp/polycore/tests/fuzz_main.cpp -o polycore_fuzz && ./polycore_fuzz 200000
// Checks: no UB / out-of-bounds on garbage input; rasterize == brute-force per-pixel
// oracle formula (in __int128) incl. int32-extreme coordinates; validate sweep == brute;
// scratch buffers return to all-zero; batch results independent of thread count.
#include <cstdio>
#include <cstdlib>
#include <random>

#include "polycore/polycore.hpp"

using namespace polycore;

static void brute_winding(const Rings& rs, i32 W, i32 H, i64 G, std::vector<i32>& w) {
  w.assign(size_t(W) * H, 0);
  for (i32 r = 0; r < H; ++r) {
    const i128 py = i128(2 * r + 1) * G;
    for (i32 ring = 0; ring < rs.nrings; ++ring) {
      const i32 b = rs.begin(ring), n = rs.size(ring);
      for (i32 k = 0; k < n; ++k) {
        const i32 v0 = b + k, v1 = b + (k + 1 == n ? 0 : k + 1);
        const i128 x0 = i128(rs.x(v0)) * 2 * W, y0 = i128(rs.y(v0)) * 2 * H;
        const i128 x1 = i128(rs.x(v1)) * 2 * W, y1 = i128(rs.y(v1)) * 2 * H;
        if (y0 == y1 || !(std::min(y0, y1) <= py && py < std::max(y0, y1))) continue;
        const int d = y1 > y0 ? 1 : -1;
        for (i32 c = 0; c < W; ++c) {
          const i128 px = i128(2 * c + 1) * G;
          const i128 num = (x0 - px) * (y1 - y0) + (py - y0) * (x1 - x0);
          if ((y1 - y0) > 0 ? num > 0 : num < 0) w[size_t(r) * W + c] += d;
        }
      }
    }
  }
}

int main(int argc, char** argv) {
  const long iters = argc > 1 ? std::atol(argv[1]) : 20000;
  std::mt19937_64 rng(12345);
  auto U = [&](i64 lo, i64 hi) { return std::uniform_int_distribution<i64>(lo, hi)(rng); };
  const i32 picks[] = {INT32_MIN, INT32_MIN + 1, -1000000000, -1, 0, 1, 499, 500, 1000, 1000000000, INT32_MAX};
  RasterScratch rsc;
  ValidateScratch vsc;
  ValidateResult v1, v2;
  std::vector<i32> wb, wo;
  long raster_checked = 0;
  for (long it = 0; it < iters; ++it) {
    const int mode = int(U(0, 3));
    const i32 R = i32(U(0, 4));
    std::vector<i32> off(R + 1, 0), xy;
    for (i32 r = 0; r < R; ++r) off[r + 1] = off[r] + i32(U(0, 9));
    for (i32 v = 0; v < 2 * off[R]; ++v) {
      i32 c;
      if (mode == 0) c = i32(U(0, 1000));
      else if (mode == 1) c = i32(U(0, 20)) * 50;
      else if (mode == 2) c = picks[U(0, 10)];
      else c = i32(U(INT32_MIN, INT32_MAX));
      xy.push_back(c);
    }
    Rings rs{xy.data(), off.data(), R};
    const i32 W = i32(U(1, 24)), H = i32(U(1, 24));
    const i64 G = U(0, 3) == 0 ? U(1, kMaxGrid) : 1000;
    std::vector<u8> eo(size_t(W) * H), nz(size_t(W) * H);
    wo.assign(size_t(W) * H, 0);
    std::vector<i32> wmin(R), wmax(R);
    RasterOut out{wo.data(), eo.data(), nz.data(), wmin.data(), wmax.data()};
    rasterize(rs, W, H, G, rsc, out);
    brute_winding(rs, W, H, G, wb);
    if (wb != wo) { std::printf("RASTER MISMATCH it=%ld\n", it); return 1; }
    for (size_t i = 0; i < wo.size(); ++i)
      if (eo[i] != u8(wo[i] & 1) || nz[i] != u8(wo[i] != 0)) { std::printf("MASK MISMATCH it=%ld\n", it); return 1; }
    // second call through the scratch winding path must give the same masks
    std::vector<u8> eo2(eo.size());
    RasterOut out2{nullptr, eo2.data(), nullptr, nullptr, nullptr};
    rasterize(rs, W, H, G, rsc, out2);
    if (eo2 != eo) { std::printf("SCRATCH PATH MISMATCH it=%ld\n", it); return 1; }
    ++raster_checked;
    // roles: random valid assignment
    std::vector<i32> roles(R, -1);
    for (i32 r = 1; r < R; ++r) if (U(0, 1)) roles[r] = 0;
    validate(rs, R ? roles.data() : nullptr, G, v1, vsc, PairMethod::Sweep);
    validate(rs, R ? roles.data() : nullptr, G, v2, vsc, PairMethod::Brute);
    if (v1.counts != v2.counts || v1.total != v2.total || v1.pair_cross != v2.pair_cross ||
        v1.pair_touch != v2.pair_touch) { std::printf("SWEEP/BRUTE MISMATCH it=%ld\n", it); return 1; }
  }
  for (i32 v : rsc.diff) if (v) { std::printf("diff scratch not zero\n"); return 1; }
  for (i32 v : rsc.wind) if (v) { std::printf("wind scratch not zero\n"); return 1; }
  // batch determinism across thread counts
  std::vector<std::vector<i32>> xys, offs;
  std::vector<Rollout> ros;
  for (int i = 0; i < 64; ++i) {
    std::vector<i32> o = {0, i32(U(3, 40))}, x;
    for (i32 v = 0; v < 2 * o[1]; ++v) x.push_back(i32(U(0, 1000)));
    xys.push_back(x); offs.push_back(o);
  }
  for (int i = 0; i < 64; ++i) {
    Rollout ro; ro.rings = Rings{xys[i].data(), offs[i].data(), 1}; ro.W = ro.H = i32(U(8, 128));
    ros.push_back(ro);
  }
  std::vector<RolloutScore> s1(64), s8(64);
  score_batch(ros.data(), 64, 1000, 1, s1.data());
  score_batch(ros.data(), 64, 1000, 8, s8.data());
  for (int i = 0; i < 64; ++i)
    if (s1[i].defects != s8[i].defects || s1[i].raster.area_nonzero != s8[i].raster.area_nonzero) {
      std::printf("BATCH MISMATCH\n"); return 1;
    }
  std::printf("ok: %ld random shapes (raster vs brute, sweep vs brute), scratch clean, batch deterministic\n", raster_checked);
  return 0;
}
